#!/usr/bin/env python3
"""Bounded full-call coloc comparisons; see README.md in this directory.

Only the fixture export command reads bundled data. Measurements accept a frozen
exchange directory, making all metadata, input arrays and component IDs explicit.
No fine-mapping is performed by these pairing measurements.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import time
from pathlib import Path

import numpy as np

import pycoloc
from pycoloc import coloc_abf, coloc_bf_bf, coloc_susie

HYPOTHESES = [f"PP.H{i}.abf" for i in range(5)]


def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return safe(value.tolist())
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, np.generic):
        return safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def export_fixture(data_dir, output, pair):
    """Export existing values bit for bit, without LD or refitting, for Python/R."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    data_dir = Path(data_dir)
    metadata = json.loads((data_dir / "metadata.json").read_text())
    fitmeta = json.loads((data_dir / "r_fits.json").read_text())
    case = {"case_id": "_".join(pair), "datasets": [], "fits": [],
            "priors": {"abf": {"p1": 1e-4, "p2": 1e-4, "p12": 1e-5},
                       "bf": {"p1": 1e-4, "p2": 1e-4, "p12": 5e-6},
                       "susie": {"p1": 1e-4, "p2": 1e-4, "p12": 5e-6}},
            "overlap_min": 0.5, "trim_by_posterior": True,
            "provenance": {"nature": "existing synthetic coloc teaching fixture",
                           "reference_lock": json.loads((data_dir / "reference.lock.json").read_text())}}

    def export_array(name, values):
        a = np.asarray(values, dtype="<f8", order="C")
        path = output / (name + ".f64")
        a.tofile(path)
        return {"file": path.name, "shape": list(a.shape),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    with np.load(data_dir / "inputs.npz", allow_pickle=False) as arrays:
        for key in pair:
            m = metadata[key]
            d = {k: m[k] for k in ("type", "N", "sdY", "snp")}
            for field in ("beta", "varbeta", "MAF"):
                d[field] = export_array(key + "_" + field, arrays[key + "_" + field])
            case["datasets"].append(d)
    with np.load(data_dir / "r_fits.npz", allow_pickle=False) as arrays:
        for key in pair:
            m = fitmeta[key]
            fit = {"variant_ids": m["snp"], "converged": m["converged"],
                   "niter": m["niter"], "null_index": None,
                   "lbf_variable": export_array(key + "_lbf_variable", arrays[key + "_lbf_variable"]),
                   "sets": {"cs": [s["members_0based"] for s in m["cs"]],
                            "cs_index": [s["component_id_0based"] for s in m["cs"]]}}
            case["fits"].append(fit)
    (output / "case.json").write_text(json.dumps(case, indent=2) + "\n")


def load_case(case, stage=None):
    """Load case.json and only stage-relevant data outside clocks.

    Array descriptors name row-major little-endian float64 files and their shape.
    Fit members and component IDs are zero based. Callers may instead supply a
    preloaded dictionary to run_case; no input preparation is hidden in run_case.
    """
    path = Path(case)
    specification = json.loads(path.read_text())

    def array(spec):
        p = path.parent / spec["file"]
        expected = int(np.prod(spec["shape"]))
        if p.stat().st_size != expected * 8:
            raise ValueError(f"Incorrect float64 byte count: {p.name}")
        return np.fromfile(p, dtype="<f8").reshape(spec["shape"])

    inputs = {k: specification[k] for k in ("case_id", "priors", "overlap_min", "trim_by_posterior")}
    if stage in (None, "abf"):
        inputs["datasets"] = [{k: array(v) if isinstance(v, dict) and "file" in v else v
                               for k, v in d.items()} for d in specification["datasets"]]
    if stage in (None, "bf", "susie"):
        inputs["fits"] = [{**fit, "lbf_variable": array(fit["lbf_variable"])}
                          for fit in specification["fits"]]
    return inputs


def run_case(inputs, stage, engine="numpy"):
    """One complete public call; returns its unmodified ColocResult."""
    priors = inputs["priors"][stage]
    if stage == "abf":
        return coloc_abf(*inputs["datasets"], **priors, engine=engine)
    f1, f2 = inputs["fits"]
    options = dict(priors, overlap_min=inputs["overlap_min"],
                   trim_by_posterior=inputs["trim_by_posterior"], engine=engine)
    if stage == "susie":
        return coloc_susie(f1, f2, **options)
    if stage == "bf":
        # Every full row, not only credible-set members or PIPs.
        return coloc_bf_bf(f1["lbf_variable"], f2["lbf_variable"],
                           snps1=f1["variant_ids"], snps2=f2["variant_ids"],
                           null_index1=f1.get("null_index"), null_index2=f2.get("null_index"),
                           **options)
    raise ValueError(f"Unknown stage: {stage}")


def snapshot(result):
    return safe({"status": result.status, "summary": result.summary.to_dict("records"),
                 "results": result.results.to_dict("records"), "priors": result.priors,
                 "metadata": result.metadata})


def compare_result(result, reference):
    """Compare raw R results without changing tolerances or hiding failed fields.

    Pair IDs are normalized from R's one-based original component indices, not
    matched by output order or PIP. Pass R's status explicitly for unavailable cases.
    The same adapter can compare native-fit pairing against official-fit pairing.
    """
    import pandas as pd

    checks = []

    def exact(field, actual, expected):
        checks.append({"field": field, "passed": safe(actual) == safe(expected),
                       "actual": safe(actual), "expected": safe(expected)})

    def numeric(field, actual, expected, rtol=0):
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        shape = a.shape == b.shape
        mask = np.isfinite(a) & np.isfinite(b) if shape else np.array([], bool)
        nonfinite = shape and np.array_equal(np.isnan(a), np.isnan(b)) and np.array_equal(np.isposinf(a), np.isposinf(b)) and np.array_equal(np.isneginf(a), np.isneginf(b))
        error = np.abs(a[mask] - b[mask]) if shape else np.array([])
        outside = int(np.sum(error > 1e-10 + rtol * np.abs(b[mask]))) if shape else None
        checks.append({"field": field, "passed": bool(shape and nonfinite and outside == 0),
                       "max_abs_error": float(error.max()) if error.size else None,
                       "n_compared": int(mask.sum()), "n_outside_tolerance": outside,
                       "shape_matches": shape, "atol": 1e-10, "rtol": rtol})

    expected_status = reference.get("status", "ok")
    exact("status", result.status, expected_status)
    if expected_status in ("no_cs", "no_overlap"):
        exact("empty summary", result.summary.empty, True)
        return {"passed": all(c["passed"] for c in checks), "checks": checks}
    rs = reference["summary"]
    if isinstance(rs, dict):
        rs = [rs]
    exact("summary row count", len(result.summary), len(rs))
    ref_table = pd.DataFrame(reference.get("results", []))
    if expected_status in ("ok", "partially_trimmed"):
        exact("reference contains SNP results", "snp" in ref_table, True)
        exact("Python contains SNP results", "snp" in result.results, True)
    if "snp" in ref_table:
        exact("SNP IDs/order", result.results.snp.tolist(), ref_table.snp.tolist())
        exact("SNP uniqueness", result.results.snp.is_unique, True)
    pair = bool(rs and "idx1" in rs[0])
    if pair:
        ref_ids = [(int(r["idx1"]) - 1, int(r["idx2"]) - 1) for r in rs]
        actual_ids = list(zip(result.summary.idx1.astype(int), result.summary.idx2.astype(int)))
        exact("original component pairs", sorted(actual_ids), sorted(ref_ids))
    for i, r in enumerate(rs):
        if pair:
            ids = ref_ids[i]
            chosen = result.summary[(result.summary.idx1 == ids[0]) & (result.summary.idx2 == ids[1])]
            if len(chosen) != 1:
                exact(f"pair {ids} occurs once", len(chosen), 1)
                continue
            ar = chosen.iloc[0]
            label = str(ids)
            refcolumn = "SNP.PP.H4.abf" if len(rs) == 1 else f"SNP.PP.H4.row{i+1}"
            actualcolumn = ar.get("snp_pp_column")
            for key in ("hit1", "hit2"):
                exact(label + " " + key, ar[key], r[key])
        else:
            if result.summary.empty:
                continue
            ar, label = result.summary.iloc[0], "ABF"
            refcolumn = actualcolumn = "SNP.PP.H4"
        exact(label + " SNP count", int(ar.nsnps), int(r["nsnps"]))
        numeric(label + " H0-H4", ar[HYPOTHESES].to_numpy(float), [r[k] for k in HYPOTHESES])
        if expected_status in ("ok", "partially_trimmed"):
            exact(label + " reference conditional column", refcolumn in ref_table, True)
            exact(label + " Python conditional column", actualcolumn in result.results, True)
        if actualcolumn in result.results and refcolumn in ref_table and "snp" in result.results:
            ref_index = ref_table.set_index("snp")
            if set(ref_index.index) == set(result.results.snp):
                numeric(label + " SNP.PP.H4", result.results[actualcolumn],
                        ref_index.loc[result.results.snp, refcolumn], rtol=1e-10)
                snp_prob = result.results[actualcolumn].to_numpy(float)
                if ar.get("h4_status", "defined") != "zero_shared_evidence":
                    exact(label + " conditional probability bounds", bool(np.isfinite(snp_prob).all() and (snp_prob >= 0).all() and (snp_prob <= 1).all()), True)
                    numeric(label + " conditional probability sum", [snp_prob.sum()], [1.0])
        if not pair and "snp" in ref_table:
            ref_index = ref_table.set_index("snp")
            if set(ref_index.index) == set(result.results.snp):
                for key in ("lABF.df1", "lABF.df2"):
                    numeric(label + " " + key, result.results[key], ref_index.loc[result.results.snp, key], rtol=1e-10)
    if result.status in ("ok", "partially_trimmed"):
        h = result.summary[HYPOTHESES].to_numpy(float)
        exact("hypothesis probability bounds", bool(np.isfinite(h).all() and (h >= 0).all() and (h <= 1).all()), True)
        numeric("hypothesis row sums", h.sum(axis=1), np.ones(len(h)))
    return {"passed": all(c["passed"] for c in checks), "checks": checks}


def environment():
    directory = Path(pycoloc.__file__).parent
    hashes = {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(directory.rglob("*.py"))}
    record = {"python": platform.python_version(), "platform": platform.platform(),
              "version": pycoloc.__version__, "source_files": hashes,
              "source_hash": hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(),
              "dependencies": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pandas")},
              "threads": {p: os.environ.get(p) for p in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")}}
    try:
        from threadpoolctl import threadpool_info
        record["blas"] = threadpool_info()
    except ImportError:
        record["blas"] = None
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--export-fixture", type=Path, help="Create a new portable input directory and exit")
    p.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "examples/data")
    p.add_argument("--pair", nargs=2, default=["D1", "D2"], choices=["D1", "D2", "D3", "D4"])
    p.add_argument("--case", type=Path)
    p.add_argument("--stage", choices=["abf", "bf", "susie"], default="abf")
    p.add_argument("--mode", choices=["validation", "timing", "memory"], default="validation")
    p.add_argument("--engine", choices=["numpy", "numba"], default="numpy")
    p.add_argument("--repeats", type=int, default=7)
    p.add_argument("--batch", type=int, default=10)
    p.add_argument("--reference", type=Path)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    if args.export_fixture:
        export_fixture(args.data_dir, args.export_fixture, args.pair)
        return 0
    if not args.case or not args.output:
        p.error("--case and --output are required for execution")
    if args.repeats < 7 or args.batch < 1:
        p.error("timing requires at least seven retained repeats and a positive batch")
    if args.output.resolve().is_relative_to(args.data_dir.resolve()):
        p.error("output must not overwrite frozen fixture data")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    record = {"case_id": json.loads(args.case.read_text())["case_id"], "stage": args.stage,
              "backend": "pycoloc", "mode": args.mode, "engine": args.engine,
              "status": "error", "records": [], "error": None}
    try:
        inputs = load_case(args.case, args.stage)
        result = None
        if args.mode == "timing":
            start = time.perf_counter()
            result = run_case(inputs, args.stage, args.engine)
            record["warmup_seconds"] = time.perf_counter() - start
            result = None
            for repeat in range(args.repeats):
                start_epoch = time.time()
                start = time.perf_counter()
                for _ in range(args.batch):
                    result = None
                    result = run_case(inputs, args.stage, args.engine)
                elapsed = time.perf_counter() - start
                end_epoch = time.time()
                record["records"].append({"repeat": repeat, "time_seconds": elapsed / args.batch,
                                          "batch_calls": args.batch, "peak_rss_mib": None,
                                          "start_epoch": start_epoch, "end_epoch": end_epoch})
        else:
            result = run_case(inputs, args.stage, args.engine)
        # Capture HWM before environment discovery or numerical serialization.
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        input_ids = [f["variant_ids"] for f in inputs["fits"]] if "fits" in inputs else [d["snp"] for d in inputs["datasets"]]
        record.update(status=result.status, n_variants=len(set(input_ids[0]) & set(input_ids[1])),
                      peak_rss_mib=peak if args.mode == "memory" else None)
        if "fits" in inputs:
            record.update(L=[len(f["lbf_variable"]) for f in inputs["fits"]],
                          niter=[f.get("niter") for f in inputs["fits"]],
                          converged=[f["converged"] for f in inputs["fits"]])
        else:
            record.update(L=None, niter=None, converged=None)
        if args.mode == "validation":
            record["result"] = snapshot(result)
            if args.reference:
                raw = json.loads(args.reference.read_text())
                record["comparison"] = compare_result(result, raw.get("result", raw))
    except Exception as exc:
        record["error"] = f"{type(exc).__name__}: {exc}"
    record["environment"] = environment()
    args.output.write_text(json.dumps(safe(record), indent=2, allow_nan=False) + "\n")
    print(json.dumps({"case_id": record["case_id"], "stage": args.stage, "mode": args.mode,
                      "status": record["status"], "error": record["error"]}))
    return int(record["error"] is not None or not record.get("comparison", {}).get("passed", True))


if __name__ == "__main__":
    raise SystemExit(main())
