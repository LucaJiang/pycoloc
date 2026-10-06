"""Stable O(m) hypothesis calculations; natural logs, float64, no fastmath."""

import numpy as np
from scipy.special import logsumexp


def _offdiag_loop(a, b):
    # Add both orientations with preceding SNPs; never subtract nearly equal
    # total and diagonal masses. This retains H3 for extreme coincident leads.
    sa = sb = total = -np.inf
    for j in range(len(a)):
        total = np.logaddexp(total, np.logaddexp(a[j] + sb, b[j] + sa))
        sa = np.logaddexp(sa, a[j])
        sb = np.logaddexp(sb, b[j])
    return total


_jit_offdiag = None


def offdiag_logsum(a, b, engine="numpy"):
    if engine == "numba":
        global _jit_offdiag
        if _jit_offdiag is None:
            try:
                from numba import njit
            except ImportError as exc:
                raise ImportError("engine='numba' requires pycoloc[jit]") from exc
            _jit_offdiag = njit(cache=False, fastmath=False)(_offdiag_loop)
        return _jit_offdiag(a, b)
    if engine != "numpy":
        raise ValueError("engine must be 'numpy' or 'numba'")
    before = np.r_[-np.inf, np.logaddexp.accumulate(b)[:-1]]
    after = np.r_[np.logaddexp.accumulate(b[::-1])[-2::-1], -np.inf]
    return logsumexp(a + np.logaddexp(before, after))


def log_weights(w, m):
    if w is None:
        return np.zeros(m)
    # Work in log space to remain invariant to enormous overall scaling.
    with np.errstate(divide="ignore"):
        lw = np.log(np.asarray(w, dtype=np.float64))
    return np.log(float(m)) + lw - logsumexp(lw)


def conditional(logmass, norm=None):
    if norm is None:
        norm = logsumexp(logmass)
    return np.exp(logmass - norm) if np.isfinite(norm) else np.full(len(logmass), np.nan)


def h4_status(log_u4, pp_h4):
    """Separate zero configuration support from a representational underflow."""
    if np.isneginf(log_u4):
        return "zero_shared_evidence"
    return "posterior_underflow" if pp_h4 == 0 else "defined"


def combine(l1, l2, priors, w1=None, w2=None, engine="numpy"):
    """Return H0–H4, conditional shared SNP PP, and log(U0), ..., log(U4).

    The shared mass uses q12 = q1*q2*p12/(p1*p2), with no further prior
    normalization. The conditional is defined by logU4, not floating PP.H4.
    """
    m = len(l1)
    lw1, lw2 = log_weights(w1, m), log_weights(w2, m)
    a = l1 + np.log(priors["p1"]) + lw1
    b = l2 + np.log(priors["p2"]) + lw2
    shared = l1 + l2 + lw1 + lw2
    shared_norm = logsumexp(shared)
    h = np.array(
        [
            0.0,
            logsumexp(a),
            logsumexp(b),
            offdiag_logsum(a, b, engine),
            np.log(priors["p12"]) + shared_norm,
        ]
    )
    return np.exp(h - logsumexp(h)), conditional(shared, shared_norm), h
