import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycoloc import coloc_bf_bf, coloc_susie
from pycoloc.result import HYPOTHESES


def test_full_rows_and_actual_component_indices(fits):
    f1, f2 = fits
    out = coloc_susie(f1, f2)
    assert out.summary[["idx1", "idx2"]].values.tolist() == [
        [0, 0],
        [2, 0],
        [0, 1],
        [2, 1],
    ]
    assert out.summary.idx1_r.tolist() == [1, 3, 1, 3]
    expected = coloc_bf_bf(
        f1["lbf_variable"][[0, 2]],
        f2["lbf_variable"],
        snps1=f1["variant_ids"],
        snps2=f2["variant_ids"],
    )
    assert_allclose(out.summary[HYPOTHESES], expected.summary[HYPOTHESES], atol=1e-15)
    assert len(out.results) == 4  # not restricted to two CS SNPs
    assert out.results.filter(regex=r"^SNP.PP.H4.row[0-9]$").shape == (4, 4)


@pytest.mark.parametrize(
    "w1,w2",
    [
        (None, None),
        ([1, 1, 1, 1], [1, 1, 1, 1]),
        ([1, 2, 3, 4], None),
        (None, [4, 3, 2, 1]),
        ([1, 2, 3, 4], [4, 3, 2, 1]),
    ],
)
def test_all_susie_weight_paths(fits, w1, w2):
    f1, f2 = fits
    out = coloc_susie(f1, f2, prior_weights1=w1, prior_weights2=w2)
    ref = coloc_bf_bf(
        f1["lbf_variable"][[0, 2]],
        f2["lbf_variable"],
        snps1=f1["variant_ids"],
        snps2=f2["variant_ids"],
        prior_weights1=w1,
        prior_weights2=w2,
    )
    assert_allclose(out.summary[HYPOTHESES], ref.summary[HYPOTHESES])


def test_reordering_and_weight_scaling(fits):
    f1, f2 = copy.deepcopy(fits)
    w1, w2 = dict(a=1, b=2, c=3, d=4), dict(a=4, b=3, c=2, d=1)
    expected = coloc_susie(f1, f2, prior_weights1=w1, prior_weights2=w2)
    order = np.array([3, 1, 0, 2])
    f2["variant_ids"] = np.asarray(f2["variant_ids"])[order]
    f2["lbf_variable"] = f2["lbf_variable"][:, order]
    f2["sets"]["cs"] = [[1], [3]]
    actual = coloc_susie(
        f1, f2, prior_weights1={k: v * 1e100 for k, v in w1.items()}, prior_weights2=w2
    )
    assert_allclose(actual.summary[HYPOTHESES], expected.summary[HYPOTHESES], atol=1e-14)


def test_null_offsets_and_full_support_overlap():
    # Three of four SNPs overlap, but virtually none of trait 1's signal does.
    a = np.array([[0, 0, 0, 20, 3.0], [2, 4, 2, 0, 1.0]])
    b = np.array([[0, 3, 1, 2.0]])
    out = coloc_bf_bf(
        a,
        b,
        snps1=["a", "b", "c", "x"],
        snps2=["a", "b", "c"],
        null_index1=4,
        null_index2=3,
    )
    assert out.status == "partially_trimmed"
    assert out.summary.idx1.tolist() == [1]
    assert out.metadata["dropped_pairs"] == [[0, 0]]
    unshift = coloc_bf_bf(
        a - a[:, -1, None],
        b - b[:, -1, None],
        snps1=["a", "b", "c", "x"],
        snps2=["a", "b", "c"],
        null_index1=4,
        null_index2=3,
    )
    assert_allclose(out.summary[HYPOTHESES], unshift.summary[HYPOTHESES])
    assert out.metadata["overlap1"][0] < 0.001
    disabled = coloc_bf_bf(
        a,
        b,
        snps1=["a", "b", "c", "x"],
        snps2=["a", "b", "c"],
        null_index1=4,
        null_index2=3,
        trim_by_posterior=False,
    )
    assert len(disabled.summary) == 2


def test_all_trimmed_null_and_no_overlap():
    with pytest.warns(UserWarning, match="posterior support"):
        out = coloc_bf_bf(
            [0, 20, 0],
            [0, 10, 0],
            snps1=["a", "x"],
            snps2=["a", "y"],
            null_index1=2,
            null_index2=2,
        )
    assert out.status == "insufficient_overlap"
    assert out.summary[HYPOTHESES[1:]].isna().all().all()
    assert out.summary["PP.H0.abf"].iloc[0] > 0
    assert out.summary.hit1.iloc[0] == "x"
    no = coloc_bf_bf([1], [2], snps1=["a"], snps2=["b"])
    assert no.status == "no_overlap" and no.summary.empty


def test_partial_overlap_weights_match_ids():
    out = coloc_bf_bf(
        [2, 4, 1],
        [0, 3, 5],
        snps1=["a", "b", "c"],
        snps2=["c", "b", "d"],
        prior_weights1=[1, 2, 3],
        prior_weights2=[3, 2, 1],
        trim_by_posterior=False,
    )
    ref = coloc_bf_bf(
        [4, 1],
        [3, 0],
        snps1=["b", "c"],
        snps2=["b", "c"],
        prior_weights1=[2, 3],
        prior_weights2=[2, 3],
    )
    assert_allclose(out.summary[HYPOTHESES], ref.summary[HYPOTHESES])


def test_empty_cs_nonconvergence_and_bad_contract(fits):
    f1, f2 = copy.deepcopy(fits)
    f1["sets"] = dict(cs=[], cs_index=[])
    with pytest.warns(UserWarning, match="no credible sets"):
        out = coloc_susie(f1, f2)
    assert out.status == "no_cs" and out.summary.empty
    f1 = copy.deepcopy(fits[0])
    f1["converged"] = False
    with pytest.warns(UserWarning, match="did not converge"):
        out = coloc_susie(f1, f2)
    assert out.metadata["converged"] == [False, True]
    f1["sets"]["cs_index"] = [0, 9]
    with pytest.raises(ValueError, match="cs_index"):
        coloc_susie(f1, f2)


@pytest.mark.parametrize(
    "bf,kw",
    [
        ([np.inf, 0], {}),
        ([np.nan, 0], {}),
        ([-np.inf, -np.inf], {}),
        ([1, 2, 3], {}),
        ([1, 2, 0], {"null_index1": 0}),
        ([1, 2, -np.inf], {"null_index1": 2}),
        ([1, 2], {"prior_weights1": [0, 0]}),
        ([1, 2], {"prior_weights1": [1, np.inf]}),
    ],
)
def test_invalid_bf_inputs(bf, kw):
    with pytest.raises(ValueError):
        coloc_bf_bf(bf, [1, 2], snps1=["a", "b"], snps2=["a", "b"], **kw)


def test_zero_bf_shared_support():
    out = coloc_bf_bf([0, -np.inf], [-np.inf, 0], snps1=["a", "b"], snps2=["a", "b"])
    assert out.summary["PP.H4.abf"].iloc[0] == 0
    assert out.results["SNP.PP.H4.abf"].isna().all()


def test_null_weight_contract_and_real_snp_named_null():
    out = coloc_bf_bf(
        [1, 2, 0],
        [2, 3, 0],
        snps1=["null", "b"],
        snps2=["null", "b"],
        null_index1=2,
        null_index2=2,
        prior_weights1={"null": 1, "b": 2},
    )
    assert out.nsnps == 2 and "null" in out.results.snp.tolist()
