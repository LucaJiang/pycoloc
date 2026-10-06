"""Independent configuration/normal-integral oracles for the 0.1.1 contract."""

import math
from decimal import Decimal, localcontext

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycoloc import coloc_abf, coloc_bf_bf, coloc_susie
from pycoloc._math import combine
from pycoloc.result import HYPOTHESES

PRIORS = dict(p1=1e-4, p2=1e-4, p12=5e-6)


def oracle(l1, l2, w1=None, w2=None, priors=PRIORS):
    """Enumerate configurations with Decimal; no production math helpers."""
    with localcontext() as ctx:
        ctx.prec = 90

        def D(x):
            return Decimal(str(x))

        m = len(l1)
        w1 = [1] * m if w1 is None else w1
        w2 = [1] * m if w2 is None else w2
        p1, p2, p12 = (D(priors[k]) for k in ("p1", "p2", "p12"))
        q1 = [m * p1 * D(w) / sum(map(D, w1)) for w in w1]
        q2 = [m * p2 * D(w) / sum(map(D, w2)) for w in w2]
        b1, b2 = [D(x).exp() for x in l1], [D(x).exp() for x in l2]
        shared = [q1[j] * q2[j] * p12 / (p1 * p2) * b1[j] * b2[j] for j in range(m)]
        u = [
            D(1),
            sum(q * b for q, b in zip(q1, b1)),
            sum(q * b for q, b in zip(q2, b2)),
            sum(
                (q1[j] * q2[k] * b1[j] * b2[k] for j in range(m) for k in range(m) if j != k), D(0)
            ),
            sum(shared),
        ]
        return (
            np.array([float(x / sum(u)) for x in u]),
            np.array([float(x / u[4]) for x in shared]) if u[4] else np.full(m, np.nan),
            np.array([float(x.ln()) if x else -np.inf for x in u]),
        )


def fit(bf, ids, null=None):
    bf = np.atleast_2d(bf)
    return dict(
        variant_ids=ids,
        lbf_variable=bf,
        converged=True,
        null_index=null,
        sets=dict(cs=[[0] for _ in bf], cs_index=list(range(len(bf)))),
    )


def public_bf(kind, l1, l2, ids1=None, ids2=None, **kwargs):
    ids1 = [str(j) for j in range(np.shape(l1)[-1])] if ids1 is None else ids1
    ids2 = ids1 if ids2 is None else ids2
    if kind == "bf":
        return coloc_bf_bf(l1, l2, snps1=ids1, snps2=ids2, **PRIORS, **kwargs)
    null1, null2 = kwargs.pop("null_index1", None), kwargs.pop("null_index2", None)
    return coloc_susie(fit(l1, ids1, null1), fit(l2, ids2, null2), **PRIORS, **kwargs)


@pytest.mark.parametrize("engine", ["numpy", "numba"])
@pytest.mark.parametrize("kind", ["abf", "bf", "susie"])
@pytest.mark.parametrize(
    "w1,w2,expected,u4",
    [
        ([9, 1], [1, 1], [0.9, 0.1], 1e-5),
        ([9, 1], [9, 1], [81 / 82, 1 / 82], 1.64e-5),
        ([1, 0], [1, 1], [1, 0], 1e-5),
        ([1, 0], [0, 1], [np.nan, np.nan], 0),
    ],
)
def test_analytic_weighted_main_all_interfaces(kind, engine, w1, w2, expected, u4):
    if kind == "abf":
        # W=V=0.04, z²=2 log(2): both normal density ratios are exactly 1.
        d = dict(
            snp=["a", "b"], type="cc", varbeta=[0.04] * 2, beta=[math.sqrt(0.08 * math.log(2))] * 2
        )
        result = coloc_abf(d, d, **PRIORS, prior_weights1=w1, prior_weights2=w2, engine=engine)
    else:
        result = public_bf(
            kind, [0.0, 0.0], [0.0, 0.0], prior_weights1=w1, prior_weights2=w2, engine=engine
        )
    row = result.summary.iloc[0]
    column = "SNP.PP.H4" if kind == "abf" else row.snp_pp_column
    assert_allclose(result.results[column], expected, atol=1e-10, rtol=1e-10)
    assert_allclose(result.results[column], result.results[column + ".weighted"], atol=0, rtol=0)
    pp, _, log_u = oracle([0.0, 0.0], [0.0, 0.0], w1, w2)
    assert_allclose(row[HYPOTHESES].to_numpy(float), pp, atol=1e-10, rtol=0)
    assert_allclose(row.logU4, math.log(u4) if u4 else -np.inf, atol=1e-10, rtol=0)
    assert_allclose(row.logU3, log_u[3], atol=1e-10, rtol=0)
    assert row.h4_status == ("defined" if u4 else "zero_shared_evidence")


@pytest.mark.parametrize("engine", ["numpy", "numba"])
@pytest.mark.parametrize("kind", ["bf", "susie"])
@pytest.mark.parametrize(
    "l1,l2,state",
    [
        ([1000.0, 0.0], [1000.0, 0.0], "defined"),
        ([2.0], [3.0], "defined"),
        ([-np.inf, 0.0], [0.0, -np.inf], "zero_shared_evidence"),
        ([-1000.0, -1001.0], [0.0, 0.0], "posterior_underflow"),
    ],
)
def test_zero_extreme_and_underflow_public(kind, engine, l1, l2, state):
    result = public_bf(kind, l1, l2, engine=engine)
    pp, snp, log_u = oracle(l1, l2)
    row = result.summary.iloc[0]
    assert_allclose(row[HYPOTHESES].to_numpy(float), pp, atol=1e-10, rtol=0)
    assert_allclose(result.results[row.snp_pp_column], snp, atol=1e-10, rtol=1e-10)
    assert_allclose([row.logU3, row.logU4], log_u[3:], atol=1e-10, rtol=0)
    assert row.h4_status == state
    if state == "posterior_underflow":
        assert row["PP.H4.abf"] == 0 and np.isfinite(row.logU4)
        assert_allclose(snp, [1 / (1 + math.exp(-1)), 1 / (1 + math.exp(1))], atol=1e-15)
    if len(l1) == 1:
        assert row["PP.H3.abf"] == 0 and np.isneginf(row.logU3)
    if l1[0] == 1000:
        assert_allclose(row.logU3, 1000 + math.log(2) + 2 * math.log(1e-4), atol=1e-10, rtol=0)


@pytest.mark.parametrize("engine", ["numpy", "numba"])
def test_abf_h4_underflow_is_defined(engine):
    v = np.nextafter(0.0, 1.0)
    d = dict(snp=["a", "b"], type="cc", beta=[0.0, 0.0], varbeta=[v, 2 * v])
    result = coloc_abf(d, d, **PRIORS, engine=engine)
    row = result.summary.iloc[0]
    assert row["PP.H4.abf"] == 0 and np.isfinite(row.logU4)
    assert row.h4_status == "posterior_underflow"
    assert_allclose(result.results["SNP.PP.H4"], [1 / 3, 2 / 3], atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("kind", ["bf", "susie"])
@pytest.mark.parametrize("engine", ["numpy", "numba"])
@pytest.mark.parametrize(
    "order1,order2", [([2, 0, 1], [0, 1, 2]), ([0, 1, 2], [1, 2, 0]), ([2, 0, 1], [1, 2, 0])]
)
@pytest.mark.parametrize("partial", [False, True])
def test_identity_symmetry_every_pair_explicit_null(kind, engine, order1, order2, partial):
    ids1 = np.array(["null", "b", "c"])
    ids2 = np.array(["null", "b", "x" if partial else "c"])
    a = np.array([[1.0, 2.0, 0.0], [0.0, 1.0, 3.0]])
    b = np.array([[2.0, 0.0, 1.0], [1.0, 3.0, 0.0]])
    w1, w2 = np.array([9.0, 1.0, 2.0]), np.array([1.0, 4.0, 3.0])

    def run(o1, o2, scale=1):
        # Shift each full component by a different explicit null BF.
        return public_bf(
            kind,
            np.c_[a[:, o1] + 4, [4.0, 4.0]],
            np.c_[b[:, o2] - 3, [-3.0, -3.0]],
            ids1[o1],
            ids2[o2],
            null_index1=3,
            null_index2=3,
            prior_weights1=w1[o1] * scale,
            prior_weights2=w2[o2],
            trim_by_posterior=False,
            engine=engine,
        )

    out = run(order1, order2, 7.0)
    base = run([0, 1, 2], [0, 1, 2])
    assert "null" in out.results.snp.to_list()
    assert_allclose(out.summary[HYPOTHESES], base.summary[HYPOTHESES], atol=1e-10, rtol=0)
    m = 2 if partial else 3
    for _, row in out.summary.iterrows():
        i, j = int(row.idx1), int(row.idx2)
        pp, snp, _ = oracle(a[i, :m], b[j, :m], w1[:m], w2[:m])
        assert_allclose(row[HYPOTHESES].to_numpy(float), pp, atol=1e-10, rtol=0)
        actual = out.results.set_index("snp").loc[ids1[:m], row.snp_pp_column]
        assert_allclose(actual, snp, atol=1e-10, rtol=1e-10)
    swapped = public_bf(
        kind,
        b,
        a,
        ids2,
        ids1,
        prior_weights1=w2,
        prior_weights2=w1,
        trim_by_posterior=False,
        engine=engine,
    )
    for _, row in base.summary.iterrows():
        reverse = swapped.summary.query("idx1 == @row.idx2 and idx2 == @row.idx1").iloc[0]
        assert_allclose(
            row[HYPOTHESES].to_numpy(float),
            reverse[[HYPOTHESES[i] for i in [0, 2, 1, 3, 4]]].to_numpy(float),
            atol=1e-10,
            rtol=0,
        )
        assert_allclose(
            base.results.set_index("snp")[row.snp_pp_column],
            swapped.results.set_index("snp").loc[base.results.snp, reverse.snp_pp_column],
            atol=1e-10,
            rtol=1e-10,
        )


@pytest.mark.parametrize("engine", ["numpy", "numba"])
@pytest.mark.parametrize(
    "order1,order2", [([2, 0, 1], [0, 1, 2]), ([0, 1, 2], [1, 2, 0]), ([2, 0, 1], [1, 2, 0])]
)
@pytest.mark.parametrize("partial", [False, True])
def test_abf_independent_permutation_and_trait_exchange(engine, order1, order2, partial):
    ids = np.array(["null", "b", "c"])
    a, b = np.array([1.0, 2.0, 0.0]), np.array([2.0, 0.0, 1.0])
    w1, w2 = np.array([9.0, 1.0, 2.0]), np.array([1.0, 4.0, 3.0])

    def dataset(logbf, order, second=False):
        names = ids.copy()
        if partial and second:
            names[2] = "x"
        return dict(
            snp=names[order],
            type="cc",
            varbeta=np.full(3, 0.04),
            beta=np.sqrt(0.08 * (math.log(2) + 2 * logbf[order])),
        )

    d1, d2 = dataset(a, order1), dataset(b, order2, True)
    out = coloc_abf(
        d1, d2, **PRIORS, prior_weights1=w1[order1] * 11, prior_weights2=w2[order2], engine=engine
    )
    m = 2 if partial else 3
    pp, snp, _ = oracle(a[:m], b[:m], w1[:m], w2[:m])
    assert_allclose(out.summary[HYPOTHESES].iloc[0], pp, atol=1e-10, rtol=0)
    assert_allclose(
        out.results.set_index("snp").loc[ids[:m], "SNP.PP.H4"], snp, atol=1e-10, rtol=1e-10
    )
    reverse = coloc_abf(
        d2, d1, **PRIORS, prior_weights1=w2[order2], prior_weights2=w1[order1], engine=engine
    )
    assert_allclose(
        reverse.summary[[HYPOTHESES[i] for i in [0, 2, 1, 3, 4]]],
        out.summary[HYPOTHESES],
        atol=1e-10,
        rtol=0,
    )
    assert_allclose(reverse.results["SNP.PP.H4"], out.results["SNP.PP.H4"], atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("kind", ["bf", "susie"])
def test_empty_intersection_weight_support_errors(kind):
    with pytest.raises(ValueError, match="zero total on the common SNP support"):
        public_bf(
            kind,
            [0.0, 0.0],
            [0.0, 0.0],
            ["a", "x"],
            ["a", "y"],
            prior_weights1=[0, 1],
            trim_by_posterior=False,
        )


@pytest.mark.parametrize("kind", ["quant", "cc"])
@pytest.mark.parametrize("sd", [0.3, 0.7])
def test_effect_prior_normal_integral_public(kind, sd):
    beta, V, sdY = 0.6, 0.04, 1.7
    W = (sd * sdY if kind == "quant" else sd) ** 2
    # Integrated N(beta;0,V+W) / null N(beta;0,V).
    expected = math.log(math.sqrt(V / (V + W)) * math.exp(beta * beta / 2 * (1 / V - 1 / (V + W))))
    d = dict(snp=["a"], type=kind, beta=[beta], varbeta=[V], sdY=sdY, pvalues=[0.5])
    out = coloc_abf(d, d, effect_priors=dict(quant=sd, cc=sd))
    assert_allclose(out.results["lABF.df1"], [expected], atol=1e-10, rtol=1e-10)
    assert_allclose(out.results["lABF.df2"], [expected], atol=1e-10, rtol=1e-10)
    assert out.metadata["dataset1"]["effect_prior_sd"] == (sd * sdY if kind == "quant" else sd)


@pytest.mark.parametrize("kind", ["quant", "cc"])
@pytest.mark.parametrize("p_side", [1, 2])
def test_explicit_effect_prior_pvalue_combination_rejected(kind, p_side):
    p = dict(snp=["a"], type=kind, pvalues=[0.05], MAF=[0.2], N=1000, s=0.3)
    beta = dict(p, beta=[0.2], varbeta=[0.04], sdY=1.7)
    args = (p, beta) if p_side == 1 else (beta, p)
    with pytest.raises(ValueError, match="effect_priors is unsupported for pvalues/MAF/N"):
        coloc_abf(*args, effect_priors=dict(quant=0.3, cc=0.4))


@pytest.mark.parametrize("engine", ["numpy", "numba"])
def test_all_negative_infinity_does_not_contaminate_hypotheses(engine):
    pp, snp, logs = combine(
        np.array([-np.inf, -np.inf]), np.array([0.0, 1.0]), PRIORS, engine=engine
    )
    expected, _, _ = oracle([-np.inf, -np.inf], [0.0, 1.0])
    assert_allclose(pp, expected, atol=1e-10, rtol=0)
    assert np.isnan(snp).all() and np.isneginf(logs[[1, 3, 4]]).all()


@pytest.mark.parametrize("kind", ["bf", "susie"])
@pytest.mark.parametrize("engine", ["numpy", "numba"])
def test_each_pair_has_its_own_zero_or_underflow_state(kind, engine):
    a = np.array([[0.0, -np.inf], [-1000.0, -1001.0]])
    b = np.array([[-np.inf, 0.0], [0.0, 0.0]])
    out = public_bf(kind, a, b, prior_weights1=[9, 1], prior_weights2=[0, 1], engine=engine)
    for _, row in out.summary.iterrows():
        pp, snp, logs = oracle(a[int(row.idx1)], b[int(row.idx2)], [9, 1], [0, 1])
        assert_allclose(row[HYPOTHESES].to_numpy(float), pp, atol=1e-10, rtol=0)
        assert_allclose(out.results[row.snp_pp_column], snp, atol=1e-10, rtol=1e-10)
        assert_allclose(out.results[row.snp_pp_column + ".weighted"], snp, atol=1e-10, rtol=1e-10)
        assert_allclose(row.logU4, logs[4], atol=1e-10, rtol=0)
        assert row.h4_status == ("zero_shared_evidence" if row.idx1 == 0 else "posterior_underflow")
    assert out.status == "ok"
