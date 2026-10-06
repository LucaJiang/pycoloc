import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.stats import norm

from pycoloc import coloc_abf
from pycoloc.abf import process_dataset
from pycoloc.result import HYPOTHESES


@pytest.mark.parametrize("kind", ["cc", "quant"])
def test_abf_formula_and_beta_precedence(datasets, kind):
    d = dict(datasets[0], type=kind, sdY=1.7)
    df, meta = process_dataset(d)
    W = (0.2 if kind == "cc" else 0.15 * 1.7) ** 2
    v = np.array(d["varbeta"])
    z = np.array(d["beta"]) / np.sqrt(v)
    assert_allclose(df.lABF, 0.5 * (np.log(v / (v + W)) + W / (v + W) * z * z), atol=1e-15)
    assert meta["branch"] == "beta/varbeta"
    other = dict(d, pvalues=[1, 1, 1, 1])
    assert_allclose(process_dataset(other)[0].lABF, df.lABF)


@pytest.mark.parametrize("kind", ["cc", "quant"])
def test_pvalue_formula(datasets, kind):
    d = dict(datasets[0], type=kind)
    del d["beta"], d["varbeta"]
    df, _ = process_dataset(d)
    f = np.array(d["MAF"])
    v = 1 / (2 * d["N"] * f * (1 - f))
    W = 0.15**2
    if kind == "cc":
        v /= d["s"] * (1 - d["s"])
        W = 0.2**2
    z = norm.isf(np.array(d["pvalues"]) / 2)
    assert_allclose(df.V, v)
    assert_allclose(df.lABF, 0.5 * (np.log(v / (v + W)) + W / (v + W) * z * z))


def test_sdy_matches_no_intercept_regression(datasets):
    d = dict(datasets[0], type="quant")
    with pytest.warns(UserWarning, match="estimating sdY"):
        _, meta = process_dataset(d)
    x = 1 / np.array(d["varbeta"])
    y = 2 * d["N"] * np.array(d["MAF"]) * (1 - np.array(d["MAF"]))
    expected = np.sqrt(np.linalg.lstsq(x[:, None], y, rcond=None)[0][0])
    assert_allclose(meta["sdY"], expected, atol=1e-14)
    assert meta["sdY_estimated"]


def test_weights_id_alignment_partial_overlap(datasets):
    d1, d2 = copy.deepcopy(datasets)
    w1, w2 = dict(a=1, b=2, c=3, d=4), dict(a=4, b=3, c=2, d=1)
    expected = coloc_abf(d1, d2, prior_weights1=w1, prior_weights2=w2)
    order = [3, 1, 0, 2]
    for d in [d1, d2]:
        for k, v in list(d.items()):
            if isinstance(v, list):
                d[k] = [v[i] for i in order]
    actual = coloc_abf(d1, d2, prior_weights1=[4, 2, 1, 3], prior_weights2=[1, 3, 4, 2])
    assert_allclose(actual.summary[HYPOTHESES], expected.summary[HYPOTHESES], atol=1e-15)
    assert actual.results.snp.tolist() == ["a", "b", "c", "d"]
    for k, v in list(d2.items()):
        if isinstance(v, list):
            d2[k] = v[:3]
    part = coloc_abf(d1, d2, prior_weights1=w1, prior_weights2={s: w2[s] for s in d2["snp"]})
    assert part.nsnps == 3
    assert part.weights["prior_weights1"].index.tolist() == ["a", "b", "d"]


def test_weighted_snp_field_corrects_reference(datasets):
    a = coloc_abf(*datasets)
    b = coloc_abf(*datasets, prior_weights1=[1, 2, 3, 4])
    expected = a.results["SNP.PP.H4"] * np.array([1, 2, 3, 4])
    assert_allclose(b.results["SNP.PP.H4"], expected / expected.sum())
    assert_allclose(b.results["SNP.PP.H4"], b.results["SNP.PP.H4.weighted"], rtol=0, atol=0)
    assert_allclose(b.results["SNP.PP.H4.weighted"].sum(), 1)
    c = coloc_abf(*datasets, prior_weights1=[0, 1, 2, 3])
    assert c.results["SNP.PP.H4.weighted"].iloc[0] == 0
    assert c.results["SNP.PP.H4"].iloc[0] == 0


@pytest.mark.parametrize(
    "bad",
    [
        [1, 2],
        [1, -1, 1, 1],
        [0, 0, 0, 0],
        [1, np.inf, 1, 1],
        [1, np.nan, 1, 1],
        {"a": 1},
    ],
)
def test_invalid_weights(datasets, bad):
    with pytest.raises(ValueError):
        coloc_abf(*datasets, prior_weights1=bad)


@pytest.mark.parametrize(
    "key,value",
    [
        ("varbeta", [0, 1, 1, 1]),
        ("varbeta", [-1, 1, 1, 1]),
        ("beta", [np.inf, 1, 1, 1]),
        ("snp", ["a", "a", "c", "d"]),
        ("MAF", [0, 0.2, 0.3, 0.4]),
        ("N", [1, 2]),
        ("s", 1),
        ("sdY", 0),
        ("pvalues", [2, 1, 1, 1]),
    ],
)
def test_invalid_statistics(datasets, key, value):
    d1, d2 = datasets
    with pytest.raises(ValueError):
        coloc_abf(dict(d1, **{key: value}), d2)


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf, 1.1, [0.1, 0.2]])
def test_invalid_base_prior(datasets, value):
    with pytest.raises(ValueError):
        coloc_abf(*datasets, p1=value)


def test_prior_adjustment_before_intersection_and_no_overlap(datasets):
    d1, d2 = datasets
    d2 = {k: (v[:2] if isinstance(v, list) else v) for k, v in d2.items()}
    with pytest.warns(UserWarning, match="nsnps"):
        out = coloc_abf(d1, d2, p1=0.4, p2=0.6, p12=0.5)
    assert out.priors == dict(p1=0.2, p2=1 / 3, p12=1 / 3)
    no = coloc_abf(d1, dict(d2, snp=["x", "y"]))
    assert no.status == "no_overlap" and no.summary.empty


def test_no_weight_mass_on_intersection(datasets):
    d1, d2 = datasets
    d2 = {k: (v[:2] if isinstance(v, list) else v) for k, v in d2.items()}
    with pytest.raises(ValueError, match="zero total"):
        coloc_abf(d1, d2, prior_weights1=[0, 0, 1, 1])


def test_effect_prior_explicit_and_stability(datasets):
    d = dict(datasets[0], type="quant", sdY=2)
    base, _ = process_dataset(d)
    changed, meta = process_dataset(d, dict(quant=0.3, cc=0.2))
    assert meta["effect_prior_sd"] == 0.6
    assert not np.allclose(base.lABF, changed.lABF)
    d = dict(snp=["a"], beta=[1e-100], varbeta=[1e-200], type="cc")
    assert np.isfinite(process_dataset(d)[0].lABF).all()
    d = dict(snp=["a"], pvalues=[np.nextafter(0.0, 1.0)], MAF=[0.2], N=1000, type="cc", s=0.2)
    assert np.isfinite(process_dataset(d)[0].lABF).all()


def test_position_conflict_and_shared_maf(datasets):
    d1, d2 = datasets
    with pytest.raises(ValueError, match="positions"):
        coloc_abf(d1, dict(d2, position=[101, 200, 300, 400]))
    a, b = copy.deepcopy(datasets)
    del a["MAF"], b["MAF"]
    assert coloc_abf(a, b, MAF=[0.1, 0.2, 0.3, 0.4]).nsnps == 4
