from decimal import Decimal, localcontext

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycoloc._math import combine, offdiag_logsum


def enumerate_hypotheses(l1, l2, p, w1=None, w2=None):
    # Independent small-fixture oracle, explicitly enumerating every model.
    with localcontext() as ctx:
        ctx.prec = 80

        def D(x):
            return Decimal(str(x))

        m = len(l1)
        w1 = [1] * m if w1 is None else w1
        w2 = [1] * m if w2 is None else w2
        q1 = [D(m) * D(p["p1"]) * D(w) / sum(map(D, w1)) for w in w1]
        q2 = [D(m) * D(p["p2"]) * D(w) / sum(map(D, w2)) for w in w2]
        b1, b2 = [D(x).exp() for x in l1], [D(x).exp() for x in l2]
        h = [
            D(1),
            sum(q * b for q, b in zip(q1, b1)),
            sum(q * b for q, b in zip(q2, b2)),
            sum(q1[i] * q2[j] * b1[i] * b2[j] for i in range(m) for j in range(m) if i != j),
            sum(
                q1[i] * q2[i] * D(p["p12"]) / (D(p["p1"]) * D(p["p2"])) * b1[i] * b2[i]
                for i in range(m)
            ),
        ]
        return np.array([float(x / sum(h)) for x in h])


@pytest.mark.parametrize(
    "w1,w2",
    [
        (None, None),
        ([1, 1, 1], [1, 1, 1]),
        ([1, 2, 4], None),
        (None, [4, 2, 1]),
        ([1, 2, 4], [4, 2, 1]),
        ([0, 1, 2], [2, 1, 0]),
    ],
)
@pytest.mark.parametrize("engine", ["numpy", "numba"])
def test_decimal_model_enumeration(w1, w2, engine):
    l1, l2 = [1.0, 4.0, -2.0], [3.0, 2.0, -1.0]
    p = dict(p1=1e-4, p2=2e-4, p12=1e-5)
    actual, _, _ = combine(np.array(l1), np.array(l2), p, w1, w2, engine)
    assert_allclose(actual, enumerate_hypotheses(l1, l2, p, w1, w2), atol=2e-15, rtol=2e-14)


def test_h3_extreme_coincident_lead_does_not_cancel():
    a, b = np.array([1000.0, 0.0]), np.array([1000.0, 0.0])
    assert_allclose(offdiag_logsum(a, b), 1000 + np.log(2), atol=1e-12)
    assert_allclose(offdiag_logsum(a, b, "numba"), 1000 + np.log(2), atol=1e-12)


def test_single_snp_h3_zero_and_disjoint_zero_weights():
    p = dict(p1=1e-4, p2=1e-4, p12=1e-5)
    pp, _, _ = combine(np.array([1000.0]), np.array([1000.0]), p)
    assert pp[3] == 0 and pp[4] == 1
    pp, weighted, _ = combine(np.array([1.0, 2.0]), np.array([2.0, 3.0]), p, [1, 0], [0, 1])
    assert pp[4] == 0 and np.isnan(weighted).all()


def test_weights_scaling_extreme_and_shared_prior_not_renormalized():
    p = dict(p1=1e-4, p2=1e-4, p12=1e-5)
    a = np.array([1.0, 2.0, 3.0])
    base = combine(a, a, p, [1, 2, 3], [3, 2, 1])[0]
    scaled = combine(a, a, p, np.array([1, 2, 3]) * 1e300, np.array([3, 2, 1]) * 1e-300)[0]
    assert_allclose(base, scaled, atol=1e-14)
    assert_allclose(base, enumerate_hypotheses(a, a, p, [1, 2, 3], [3, 2, 1]), atol=1e-14)


def test_invalid_engine():
    with pytest.raises(ValueError, match="engine"):
        offdiag_logsum(np.array([1.0]), np.array([1.0]), "fast")
