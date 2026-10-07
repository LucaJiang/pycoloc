"""Integrity and false-pass regression checks for the comparison adapter."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("coloc_case", ROOT / "benchmarks/coloc_case.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def test_export_preserves_bits_and_loads_only_stage_inputs(tmp_path):
    path = tmp_path / "exchange"
    runner.export_fixture(ROOT / "examples/data", path, ["D1", "D2"])
    abf = runner.load_case(path / "case.json", "abf")
    pairing = runner.load_case(path / "case.json", "susie")
    assert "fits" not in abf and "datasets" not in pairing
    assert all("LD" not in d for d in abf["datasets"])
    with np.load(ROOT / "examples/data/inputs.npz", allow_pickle=False) as original:
        for key, dataset in zip(["D1", "D2"], abf["datasets"]):
            for field in ["beta", "varbeta", "MAF"]:
                assert np.array_equal(dataset[field], original[key + "_" + field])
    with np.load(ROOT / "examples/data/r_fits.npz", allow_pickle=False) as original:
        for key, fit in zip(["D1", "D2"], pairing["fits"]):
            assert np.array_equal(fit["lbf_variable"], original[key + "_lbf_variable"])


def test_comparison_rejects_wrong_order_probabilities_missing_columns_and_pairs(tmp_path):
    path = tmp_path / "exchange"
    runner.export_fixture(ROOT / "examples/data", path, ["D3", "D4"])
    result = runner.run_case(runner.load_case(path / "case.json", "susie"), "susie")
    reference = json.loads((ROOT / "examples/data/r_coloc.json").read_text())["cases"]["susie_D3_D4"]["result"]
    assert runner.compare_result(result, reference)["passed"]
    changed = deepcopy(reference)
    changed["results"] = changed["results"][::-1]
    assert not runner.compare_result(result, changed)["passed"]
    changed = deepcopy(reference)
    changed["summary"][0]["PP.H4.abf"] += 1e-4
    assert not runner.compare_result(result, changed)["passed"]
    changed = deepcopy(reference)
    changed["summary"] = changed["summary"][:1]
    assert not runner.compare_result(result, changed)["passed"]
    changed = deepcopy(reference)
    for row in changed["results"]:
        row.pop("SNP.PP.H4.row1")
    assert not runner.compare_result(result, changed)["passed"]


def test_unavailable_states_remain_distinct():
    result = runner.coloc_bf_bf([0.0], [0.0], snps1=["a"], snps2=["b"])
    assert runner.compare_result(result, {"status": "no_overlap"})["passed"]
    assert not runner.compare_result(result, {"status": "no_cs"})["passed"]
