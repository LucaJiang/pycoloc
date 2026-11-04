"""Consume actual optional prusie results, including unavailable and diagnostic fits."""

import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycoloc import coloc_bf_bf, coloc_susie
from pycoloc.result import HYPOTHESES

prusie = pytest.importorskip("prusie")
DATA = Path(__file__).resolve().parents[1] / "examples/data"


def fit(dataset, **overrides):
    meta = json.loads((DATA / "metadata.json").read_text())[dataset]
    parameters = json.loads((DATA / "parameters.json").read_text())[dataset]
    params = {
        k: v
        for k, v in parameters.items()
        if k not in {"dataset", "input_kind", "n", "z_definition"}
    }
    params.update(overrides)
    with np.load(DATA / "inputs.npz", allow_pickle=False) as arrays:
        return prusie.susie_rss(
            z=arrays[dataset + "_beta"] / np.sqrt(arrays[dataset + "_varbeta"]),
            R=arrays[dataset + "_LD"],
            n=meta["N"],
            variant_ids=meta["snp"],
            **params,
        )


@pytest.mark.parametrize("refine", [False, True])
def test_native_results_preserve_full_rows_and_component_ids(refine):
    left, right = fit("D3", refine=refine), fit("D4", refine=refine)
    assert left.converged and right.converged
    actual = coloc_susie(left, right)
    ids1, ids2 = np.asarray(left.sets["cs_index"]), np.asarray(right.sets["cs_index"])
    expected = coloc_bf_bf(
        left.lbf_variable[ids1],
        right.lbf_variable[ids2],
        snps1=left.variant_ids,
        snps2=right.variant_ids,
    )
    assert actual.status == "ok"
    assert len(actual.results) == 500
    assert actual.summary.idx1.tolist() == ids1[expected.summary.idx1].tolist()
    assert actual.summary.idx2.tolist() == ids2[expected.summary.idx2].tolist()
    assert_allclose(actual.summary[HYPOTHESES], expected.summary[HYPOTHESES], atol=0, rtol=0)
    assert actual.metadata["converged"] == [True, True]


def test_native_no_cs_is_unavailable():
    left, right = fit("D3", min_abs_corr=1), fit("D4")
    assert len(left.sets["cs_index"]) == 0
    with pytest.warns(UserWarning, match="no credible sets"):
        actual = coloc_susie(left, right)
    assert actual.status == "no_cs"
    assert actual.summary.empty and actual.results.empty


def test_native_nonconvergence_stays_visible():
    left, right = fit("D3", max_iter=1), fit("D4")
    assert not left.converged
    with pytest.warns(UserWarning, match="did not converge"):
        actual = coloc_susie(left, right)
    assert actual.metadata["converged"] == [False, True]
