#!/usr/bin/env python3
"""Offline public example: check real pycoloc calls against frozen R and mathematical expectations.

Run from a checkout/source archive after installing pycoloc. No R or prusie is
needed unless --native-pair is explicitly requested. No expectation is regenerated.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import warnings
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np
import pandas as pd

import pycoloc
from pycoloc import coloc_abf, coloc_bf_bf, coloc_susie

HYPOTHESES = [f"PP.H{i}.abf" for i in range(5)]
PRIORS = dict(p1=1e-4, p2=1e-4, p12=5e-6)
PAIRS = (("D1", "D2"), ("D3", "D4"))


def read_json(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, np.ndarray) and value.ndim == 0:
        return json_safe(value.item())
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


class Checks:
    def __init__(self):
        self.rows = []
        self.diagnostics = []

    def exact(self, case, field, actual, expected):
        ok = json_safe(actual) == json_safe(expected)
        self.rows.append(
            dict(
                case=case,
                field=field,
                passed=bool(ok),
                kind="structural",
                actual=json_safe(actual),
                expected=json_safe(expected),
            )
        )
        return ok

    def numeric(
        self,
        case,
        field,
        actual,
        expected,
        *,
        atol=1e-10,
        rtol=0,
        reference="official R",
        diagnostic=False,
    ):
        a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        shape = a.shape == b.shape
        finite = (
            shape
            and np.array_equal(np.isfinite(a), np.isfinite(b))
            and np.array_equal(np.isposinf(a), np.isposinf(b))
            and np.array_equal(np.isneginf(a), np.isneginf(b))
        )
        missing = shape and np.array_equal(np.isnan(a), np.isnan(b))
        if shape:
            mask = np.isfinite(a) & np.isfinite(b)
            errors = np.abs(a[mask] - b[mask])
            bound = atol + rtol * np.abs(b[mask])
            ok = finite and missing and bool(np.all(errors <= bound))
            maximum = float(errors.max()) if errors.size else None
            outside = int(np.sum(errors > bound))
        else:
            ok, maximum, outside = False, None, None
        row = dict(
            case=case,
            field=field,
            passed=bool(ok),
            kind="numeric",
            reference=reference,
            max_abs_error=maximum,
            n_compared=int(a.size),
            n_outside_tolerance=outside,
            shape_matches=shape,
            finite_mask_matches=bool(finite),
            missing_mask_matches=bool(missing),
            atol=atol,
            rtol=rtol,
        )
        (self.diagnostics if diagnostic else self.rows).append(row)
        return row

    def probability(self, case, field, values, *, rows=False):
        a = np.asarray(values, float)
        self.exact(
            case,
            field + " finite in [0,1]",
            bool(np.all(np.isfinite(a) & (a >= 0) & (a <= 1))),
            True,
        )
        self.numeric(
            case,
            field + " normalized",
            a.sum(axis=1) if rows else [a.sum()],
            np.ones(len(a)) if rows else [1.0],
            reference="probability validity",
        )


def verify_integrity(data_dir):
    manifest = read_json(data_dir / "manifest.json")
    files = read_json(data_dir / "checksums.json")["files"]
    native_files = read_json(data_dir / "native_checksums.json")["files"]
    if set(native_files) != {"expected_native.json"}:
        raise ValueError("Native checksum record must identify exactly expected_native.json")
    default_files = read_json(data_dir / "default_abf/checksums.json")["files"]
    files = {
        **files,
        **native_files,
        **{"default_abf/" + name: entry for name, entry in default_files.items()},
    }
    required = {
        "inputs.npz",
        "metadata.json",
        "r_fits.npz",
        "r_fits.json",
        "r_coloc.json",
        "parameters.json",
        "reference.lock.json",
        "independent_oracles.json",
        "expected_native.json",
        "variants.tsv",
        "default_abf/r_default_abf.json",
    }
    if not required.issubset(files):
        raise ValueError(
            "Manifest is missing required input/reference checksums: "
            + ", ".join(sorted(required - set(files)))
        )
    if not files:
        raise ValueError("The data manifest contains no checksums")
    count = 0
    for name, entry in files.items():
        expected = entry["sha256"] if isinstance(entry, dict) else entry
        path = (data_dir / name).resolve()
        if not path.is_relative_to(data_dir.resolve()) or not path.is_file():
            raise ValueError(f"Missing or unsafe manifest file: {name}")
        actual = sha256(path)
        if actual != expected:
            raise ValueError(
                f"Checksum mismatch for {name}: expected {expected}, got {actual}. Restore the frozen bundle; goldens are never regenerated."
            )
        count += 1
    return manifest, count


def load_data(data_dir):
    metadata = read_json(data_dir / "metadata.json")
    datasets = {}
    with np.load(data_dir / "inputs.npz", allow_pickle=False) as arrays:
        for key in ("D1", "D2", "D3", "D4"):
            m = metadata[key]
            datasets[key] = dict(
                snp=np.asarray(m["snp"], dtype=str),
                type=m["type"],
                N=m["N"],
                sdY=m["sdY"],
                beta=arrays[key + "_beta"],
                varbeta=arrays[key + "_varbeta"],
                MAF=arrays[key + "_MAF"],
                LD=arrays[key + "_LD"],
            )
    fits = {}
    fitmeta = read_json(data_dir / "r_fits.json")
    with np.load(data_dir / "r_fits.npz", allow_pickle=False) as arrays:
        for key, meta in fitmeta.items():
            fits[key] = dict(
                variant_ids=meta["snp"],
                lbf_variable=arrays[key + "_lbf_variable"],
                alpha=arrays[key + "_alpha"],
                pip=arrays[key + "_pip"],
                converged=meta["converged"],
                null_index=None,
                sets=dict(
                    cs=[s["members_0based"] for s in meta["cs"]],
                    cs_index=[s["component_id_0based"] for s in meta["cs"]],
                ),
            )
    return datasets, fits, fitmeta


def decimal_oracle(logbf1, logbf2, w1=None, w2=None):
    """Independent high-precision configuration sums; no pycoloc math helper.

    The two-SNP counterexample enumerates all distinct configurations. Larger
    fixed fixtures use positive off-diagonal row sums at 80-digit precision.
    """
    with localcontext() as ctx:
        ctx.prec = 80

        def D(x):
            return Decimal(str(x))

        m = len(logbf1)
        w1 = [1] * m if w1 is None else w1
        w2 = [1] * m if w2 is None else w2
        p1, p2, p12 = (D(PRIORS[k]) for k in ("p1", "p2", "p12"))
        q1 = [m * p1 * D(x) / sum(map(D, w1)) for x in w1]
        q2 = [m * p2 * D(x) / sum(map(D, w2)) for x in w2]
        a = [D(x).exp() for x in logbf1]
        b = [D(x).exp() for x in logbf2]
        a1, a2 = [x * y for x, y in zip(q1, a)], [x * y for x, y in zip(q2, b)]
        shared = [q1[j] * q2[j] * p12 / (p1 * p2) * a[j] * b[j] for j in range(m)]
        if m <= 4:
            offdiag = sum((a1[j] * a2[k] for j in range(m) for k in range(m) if j != k), D(0))
        else:
            total2 = sum(a2)
            offdiag = sum(x * (total2 - a2[j]) for j, x in enumerate(a1))
        u = [D(1), sum(a1), sum(a2), offdiag, sum(shared)]
        h = np.array([float(x / sum(u)) for x in u])
        snp = np.array([float(x / u[4]) for x in shared]) if u[4] else np.full(m, np.nan)
        return h, snp


def snapshot(result):
    return json_safe(
        dict(
            status=result.status,
            summary=result.summary.to_dict("records"),
            results=result.results.to_dict("records"),
        )
    )


def compare_abf(checks, key, result, reference, *, diagnostic=False):
    summary = np.asarray([reference["summary"][k] for k in ["nsnps", *HYPOTHESES]])
    if not diagnostic:
        checks.exact(key, "SNP count", result.nsnps, int(summary[0]))
    checks.numeric(
        key, "H0–H4", result.summary[HYPOTHESES].to_numpy()[0], summary[1:], diagnostic=diagnostic
    )
    ref = pd.DataFrame(reference["results"]).set_index("snp")
    actual = result.results.set_index("snp")
    if not diagnostic:
        checks.exact(key, "SNP identities/order", actual.index.to_list(), ref.index.to_list())
    ref = ref.loc[actual.index]
    for column in ("lABF.df1", "lABF.df2", "SNP.PP.H4"):
        checks.numeric(key, column, actual[column], ref[column], rtol=1e-10, diagnostic=diagnostic)
    if not diagnostic:
        checks.probability(key, "H0–H4", result.summary[HYPOTHESES], rows=True)
        checks.probability(key, "SNP.PP.H4", actual["SNP.PP.H4"])


def compare_pair(checks, key, result, reference, *, diagnostic=False):
    rows = reference["summary"]
    expected_ids = [(int(r["idx1"]) - 1, int(r["idx2"]) - 1) for r in rows]
    actual_ids = [(int(r.idx1), int(r.idx2)) for r in result.summary.itertuples()]
    if not diagnostic:
        checks.exact(key, "original component IDs", sorted(actual_ids), sorted(expected_ids))
        checks.exact(key, "result status", result.status, "ok")
    refsnps = pd.DataFrame(reference["results"]).set_index("snp")
    actualsnps = result.results.set_index("snp")
    if not diagnostic:
        checks.exact(
            key, "SNP identities/order", actualsnps.index.to_list(), refsnps.index.to_list()
        )
    for i, (pair, row) in enumerate(zip(expected_ids, rows)):
        selected = result.summary.query("idx1 == @pair[0] and idx2 == @pair[1]")
        if selected.empty:
            checks.diagnostics.append(dict(case=key, field="missing reference pair", pair=pair))
            continue
        ar = selected.iloc[0]
        refcolumn = "SNP.PP.H4.abf" if len(rows) == 1 else f"SNP.PP.H4.row{i + 1}"
        label = f"component {pair[0]}/{pair[1]}"
        checks.numeric(
            key,
            label + " H0–H4",
            ar[HYPOTHESES].to_numpy(float),
            [row[k] for k in HYPOTHESES],
            diagnostic=diagnostic,
        )
        checks.numeric(
            key,
            label + " SNP.PP.H4",
            actualsnps[ar.snp_pp_column],
            refsnps.loc[actualsnps.index, refcolumn],
            rtol=1e-10,
            diagnostic=diagnostic,
        )
        if not diagnostic:
            checks.exact(key, label + " lead SNPs", [ar.hit1, ar.hit2], [row["hit1"], row["hit2"]])
            checks.exact(key, label + " SNP count", int(ar.nsnps), int(row["nsnps"]))
            checks.probability(key, label + " hypotheses", ar[HYPOTHESES].to_numpy(float))
            checks.probability(key, label + " SNP posterior", actualsnps[ar.snp_pp_column])


def evaluate(data_dir, *, native_pair=False, expected_native=None):
    checks = Checks()
    datasets, fits, fitmeta = load_data(data_dir)
    refs = read_json(data_dir / "r_coloc.json")
    default_refs = read_json(data_dir / "default_abf/r_default_abf.json")
    parameters = read_json(data_dir / "parameters.json")
    frozen_oracles = read_json(data_dir / "independent_oracles.json")["cases"]
    actual_results = {}
    order_manifest = pd.read_csv(data_dir / "variants.tsv", sep="\t")
    for key, d in datasets.items():
        checks.exact(key, "500 unique variants", len(d["snp"]) == len(set(d["snp"])) == 500, True)
        checks.exact(
            key,
            "finite signed effects; positive variance",
            bool(
                np.isfinite(d["beta"]).all()
                and np.isfinite(d["varbeta"]).all()
                and (d["varbeta"] > 0).all()
            ),
            True,
        )
        declared = order_manifest.loc[order_manifest.dataset == key]
        checks.exact(
            key, "variant manifest identities", declared.variant_id.tolist(), d["snp"].tolist()
        )
        for column in ("array_index_0based", "ld_row_0based", "ld_column_0based"):
            checks.exact(
                key, "variant manifest " + column, declared[column].tolist(), list(range(500))
            )
        checks.exact(
            key,
            "float64 input arrays",
            all(d[f].dtype == np.float64 for f in ("beta", "varbeta", "LD", "MAF")),
            True,
        )
        checks.exact(
            key,
            "finite bounded signed LD",
            bool(np.isfinite(d["LD"]).all() and (np.abs(d["LD"]) <= 1).all()),
            True,
        )
        checks.exact(key, "LD shape/order", d["LD"].shape, (len(d["snp"]), len(d["snp"])))
        checks.numeric(
            key, "LD symmetry", d["LD"], d["LD"].T, atol=1e-14, reference="input validity"
        )
        checks.numeric(
            key,
            "LD unit diagonal",
            np.diag(d["LD"]),
            np.ones(500),
            atol=1e-14,
            reference="input validity",
        )
        checks.exact(key, "fit SNP order", fits[key]["variant_ids"], d["snp"].tolist())
        checks.probability(key, "frozen fit alpha rows", fits[key]["alpha"], rows=True)
        for cs in fitmeta[key]["cs"]:
            actual_coverage = fits[key]["alpha"][
                cs["component_id_0based"], cs["members_0based"]
            ].sum()
            checks.numeric(
                key,
                "actual returned coverage component " + str(cs["component_id_0based"]),
                [actual_coverage],
                [cs["coverage"]],
                reference="actual R CS coverage / alpha mass",
            )
    for left, right in PAIRS:
        name = left + "_" + right
        default_key = "abf_default_" + name
        default_result = coloc_abf(datasets[left], datasets[right])
        compare_abf(
            checks, default_key, default_result, default_refs["cases"][default_key]["result"]
        )
        checks.exact(
            default_key,
            "literal public default priors",
            default_result.priors,
            default_refs["cases"][default_key]["actual_default_priors"],
        )
        actual_results[default_key] = snapshot(default_result)
        abf = coloc_abf(datasets[left], datasets[right], **PRIORS)
        compare_abf(checks, "abf_" + name, abf, refs["cases"]["abf_" + name]["result"])
        actual_results["abf_" + name] = snapshot(abf)
        result = coloc_susie(
            fits[left], fits[right], **PRIORS, overlap_min=0.5, trim_by_posterior=True
        )
        compare_pair(checks, "susie_" + name, result, refs["cases"]["susie_" + name]["result"])
        actual_results["susie_" + name] = snapshot(result)
        # A second public API consumes precisely the same complete component BF rows.
        bf = coloc_bf_bf(
            fits[left]["lbf_variable"][fits[left]["sets"]["cs_index"]],
            fits[right]["lbf_variable"][fits[right]["sets"]["cs_index"]],
            snps1=fits[left]["variant_ids"],
            snps2=fits[right]["variant_ids"],
            **PRIORS,
        )
        checks.numeric(
            "bf_pairing_" + name,
            "same full BF rows H0–H4",
            bf.summary[HYPOTHESES],
            result.summary[HYPOTHESES],
            reference="public API consistency (shared production implementation)",
        )
    # No-CS is an unavailable calculation, not a numerical zero for H4.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        no_cs = coloc_susie(fits["D3_strict_purity"], fits["D4"], **PRIORS)
    checks.exact("no_cs", "status", no_cs.status, "no_cs")
    checks.exact("no_cs", "no fabricated posterior rows", len(no_cs.summary), 0)
    checks.exact(
        "no_cs",
        "reference unavailable",
        refs["cases"]["susie_strict_purity_D4"]["result"],
        [{"nsnps": None}],
    )
    checks.exact(
        "no_cs",
        "actual R coverage unavailable",
        fitmeta["D3_strict_purity"]["cs_coverage_available"],
        False,
    )
    actual_results["no_cs"] = snapshot(no_cs)
    # Independently enumerate the small fixed weighted counterexample.
    for name, w1, w2 in (
        ("weighted_counterexample", [9, 1], [1, 1]),
        ("partial_zero", [1, 0], [1, 1]),
        ("disjoint_support", [1, 0], [0, 1]),
    ):
        r = coloc_bf_bf(
            [0.0, 0.0],
            [0.0, 0.0],
            snps1=["a", "b"],
            snps2=["a", "b"],
            prior_weights1=w1,
            prior_weights2=w2,
            **PRIORS,
        )
        h, snp = decimal_oracle([0.0, 0.0], [0.0, 0.0], w1, w2)
        checks.numeric(
            name,
            "H0–H4",
            r.summary[HYPOTHESES].iloc[0],
            h,
            reference="80-digit Decimal configuration oracle",
        )
        checks.numeric(
            name,
            "SNP.PP.H4",
            r.results["SNP.PP.H4.abf"],
            snp,
            reference="80-digit Decimal configuration oracle",
        )
        checks.exact(
            name,
            "H4 state",
            r.summary.h4_status.iloc[0],
            "zero_shared_evidence" if name == "disjoint_support" else "defined",
        )
        if name in ("weighted_counterexample", "partial_zero"):
            oracle_key = (
                "weighted_equal_bf_counterexample"
                if name == "weighted_counterexample"
                else "partial_zero_bf_counterexample"
            )
            frozen = frozen_oracles[oracle_key]
            checks.numeric(
                name,
                "frozen independent hypotheses",
                r.summary[HYPOTHESES].iloc[0],
                frozen["hypothesis_posterior"],
                reference="independent frozen Decimal oracle",
            )
            checks.numeric(
                name,
                "frozen independent SNP posterior",
                r.results["SNP.PP.H4.abf"],
                frozen["conditional_shared_posterior"],
                reference="independent frozen Decimal oracle",
            )
        actual_results[name] = snapshot(r)
        if name == "weighted_counterexample":
            raw = refs["cases"]["weighted_equal_bf_counterexample"]["result"]
            compare_pair(checks, name + "_raw_R", r, raw, diagnostic=True)
    checks.exact(
        "partial_zero",
        "raw R BF wrapper rejects zeros",
        bool(refs["cases"]["partial_zero_bf_counterexample"]["error"]),
        True,
    )
    # Fixed, outcome-independent 500-SNP weights expose the upstream sort/weight mismatch.
    w1, w2 = 1 + np.arange(1, 501) % 5, 1 + np.arange(1, 501) % 3
    with warnings.catch_warnings(record=True) as weight_warnings:
        warnings.simplefilter("always")
        weighted = coloc_abf(
            datasets["D1"], datasets["D2"], prior_weights1=w1, prior_weights2=w2, **PRIORS
        )
    tab = weighted.results.set_index("snp")
    ids = datasets["D1"]["snp"].tolist()
    h, snp = decimal_oracle(tab.loc[ids, "lABF.df1"], tab.loc[ids, "lABF.df2"], w1, w2)
    checks.numeric(
        "abf_weighted_500",
        "H0–H4",
        weighted.summary[HYPOTHESES].iloc[0],
        h,
        reference="80-digit Decimal configuration oracle",
    )
    checks.numeric(
        "abf_weighted_500",
        "SNP.PP.H4",
        tab.loc[ids, "SNP.PP.H4"],
        snp,
        reference="80-digit Decimal configuration oracle",
    )
    compare_abf(
        checks,
        "abf_weighted_500_raw_R",
        weighted,
        refs["cases"]["abf_D1_D2_weighted"]["result"],
        diagnostic=True,
    )
    actual_results["abf_weighted_500"] = snapshot(weighted)
    # Demonstrate propagation of the actual, explicitly nonconverged frozen fit.
    with warnings.catch_warnings(record=True) as nonconvergence_warnings:
        warnings.simplefilter("always")
        unconverged = coloc_susie(fits["D3_one_iteration"], fits["D4"], **PRIORS)
    checks.exact(
        "one_iteration", "nonconvergence retained", unconverged.metadata["converged"], [False, True]
    )
    checks.exact(
        "one_iteration",
        "nonconvergence warning emitted",
        any("did not converge" in str(w.message) for w in nonconvergence_warnings),
        True,
    )
    # Frozen native results are a separate regression checkpoint, never R goldens.
    if expected_native is not None:
        for name, expected in expected_native["cases"].items():
            current = actual_results[name]
            checks.exact(name + "_native_golden", "status", current["status"], expected["status"])
            a, b = pd.DataFrame(current["summary"]), pd.DataFrame(expected["summary"])
            checks.exact(name + "_native_golden", "row count", len(a), len(b))
            if len(a) and len(b):
                checks.numeric(
                    name + "_native_golden",
                    "H0–H4",
                    a[HYPOTHESES],
                    b[HYPOTHESES],
                    reference="frozen pycoloc checkpoint",
                )
            a, b = pd.DataFrame(current["results"]), pd.DataFrame(expected["results"])
            if len(a) and len(b):
                checks.exact(name + "_native_golden", "SNP IDs", a.snp.to_list(), b.snp.to_list())
                for col in b.columns:
                    if col.startswith("SNP.PP.H4"):
                        checks.numeric(
                            name + "_native_golden",
                            col,
                            a[col],
                            b[col],
                            rtol=1e-10,
                            reference="frozen pycoloc checkpoint",
                        )
    native_version = None
    native_fit_details = {}
    if native_pair:
        import prusie  # Optional, explicit pair mode only.

        native_version = prusie.__version__
        native_fits = {}
        for key in ("D1", "D2", "D3", "D4"):
            d = datasets[key]
            p = {
                k: v
                for k, v in parameters[key].items()
                if k not in ("dataset", "input_kind", "n", "z_definition")
            }
            native_fits[key] = prusie.susie_rss(
                z=d["beta"] / np.sqrt(d["varbeta"]), R=d["LD"], n=d["N"], variant_ids=d["snp"], **p
            )
            f = native_fits[key]
            native_fit_details[key] = json_safe(
                dict(
                    converged=f.converged,
                    niter=f.niter,
                    sets=f.sets,
                    parameters=f.params,
                    backend_version=f.backend_version,
                    warnings=f.warnings,
                )
            )
            for field in ("alpha", "lbf_variable"):
                checks.numeric(
                    "native_" + key,
                    field + " raw diagnostic",
                    getattr(f, field),
                    fits[key][field],
                    reference="official susieR 0.16.6",
                    diagnostic=True,
                )
            checks.numeric(
                "native_" + key,
                "PIP",
                native_fits[key].pip,
                fits[key]["pip"],
                atol=1e-5,
                reference="official susieR 0.16.6",
            )
            checks.exact(
                "native_" + key,
                "PIP probability validity",
                bool(
                    np.isfinite(native_fits[key].pip).all()
                    and (native_fits[key].pip >= 0).all()
                    and (native_fits[key].pip <= 1).all()
                ),
                True,
            )
            checks.probability("native_" + key, "alpha rows", native_fits[key].alpha, rows=True)
            checks.exact(
                "native_" + key,
                "positive residual / nonnegative prior variances",
                bool(
                    np.isfinite(native_fits[key].sigma2)
                    and native_fits[key].sigma2 > 0
                    and np.isfinite(native_fits[key].V).all()
                    and (np.asarray(native_fits[key].V) >= 0).all()
                ),
                True,
            )
            checks.diagnostics.append(
                dict(
                    case="native_" + key,
                    field="convergence",
                    actual=bool(native_fits[key].converged),
                    reference=fitmeta[key]["converged"],
                    niter=int(native_fits[key].niter),
                    reference_niter=fitmeta[key]["niter"],
                )
            )
        for left, right in PAIRS:
            name = "susie_" + left + "_" + right
            r = coloc_susie(native_fits[left], native_fits[right], **PRIORS)
            checks.exact("native_" + name, "computed status", r.status, "ok")
            checks.probability("native_" + name, "H0–H4 validity", r.summary[HYPOTHESES], rows=True)
            compare_pair(
                checks, "native_" + name, r, refs["cases"][name]["result"], diagnostic=True
            )
            actual_results["native_" + name] = snapshot(r)
    return dict(
        passed=all(row["passed"] for row in checks.rows),
        checks=checks.rows,
        diagnostics=checks.diagnostics,
        actual_results=actual_results,
        warnings=[str(w.message) for w in caught]
        + [str(w.message) for w in nonconvergence_warnings]
        + [str(w.message) for w in weight_warnings],
        native_pair=native_pair,
        prusie_version=native_version,
        native_fit_details=native_fit_details,
        native_pair_gate="PIP atol 1e-5/rtol 0 plus validity; downstream differences are diagnostics",
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Write report.json and actual_results.json here (never inside frozen data)",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
        help="Frozen bundle location; default: examples/data next to this script",
    )
    parser.add_argument(
        "--native-pair",
        action="store_true",
        help="Additionally fit with separately installed prusie; default standalone mode never imports it",
    )
    args = parser.parse_args(argv)
    data, output = args.data_dir.resolve(), args.output_dir.resolve()
    if output.is_relative_to(data):
        parser.error("--output-dir must be outside the frozen data directory")
    output.mkdir(parents=True, exist_ok=True)
    report = dict(
        passed=False,
        environment=dict(
            python=platform.python_version(),
            platform=platform.platform(),
            pycoloc=pycoloc.__version__,
            backend="NumPy/SciPy",
            dependencies={
                p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pandas", "matplotlib")
            },
        ),
        command=["python", "examples/check_example.py", "--output-dir", str(args.output_dir)]
        + (["--native-pair"] if args.native_pair else []),
    )
    try:
        manifest, count = verify_integrity(data)
        expected = read_json(data / "expected_native.json")
        report.update(
            manifest_sha256=sha256(data / "manifest.json"),
            reference_lock=read_json(data / "reference.lock.json"),
            runner_sha256=sha256(Path(__file__)),
            integrity_files_checked=count,
            fixture_snp_count=manifest["snp_count"],
            expected_pycoloc_version=expected["pycoloc_version"],
            version_differs=expected["pycoloc_version"] != pycoloc.__version__,
        )
        report.update(evaluate(data, native_pair=args.native_pair, expected_native=expected))
    except Exception as exc:
        report.update(passed=False, error=f"{type(exc).__name__}: {exc}")
    actual = report.pop("actual_results", {})
    (output / "actual_results.json").write_text(
        json.dumps(json_safe(actual), indent=2, allow_nan=False) + "\n"
    )
    (output / "report.json").write_text(
        json.dumps(json_safe(report), indent=2, allow_nan=False) + "\n"
    )
    print(f"pycoloc {pycoloc.__version__} · Python {platform.python_version()} · NumPy/SciPy")
    if report.get("version_differs"):
        print(
            f"NOTE: frozen checkpoint used pycoloc {report['expected_pycoloc_version']}; version difference alone is not a failure."
        )
    for row in report.get("checks", []):
        if row.get("kind") == "numeric":
            error = row["max_abs_error"]
            print(
                f"{'PASS' if row['passed'] else 'FAIL'} {row['case']} / {row['field']}: max |Δ|={error}; atol={row['atol']}, rtol={row['rtol']}"
            )
        elif not row["passed"]:
            print(
                f"FAIL {row['case']} / {row['field']}: actual={row['actual']}, expected={row['expected']}"
            )
    for message in report.get("warnings", []):
        print("NOTE: " + message)
    if "error" in report:
        print("FAIL " + report["error"])
    print(
        f"{'PASS' if report['passed'] else 'FAIL'}: {len(report.get('checks', []))} checks; {len(report.get('diagnostics', []))} separately recorded diagnostics. Report: {output / 'report.json'}"
    )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
