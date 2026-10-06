"""Pair full component Bayes-factor rows, retaining original component IDs."""

import warnings
from collections.abc import Mapping

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from . import _validation as val
from ._math import combine, h4_status
from .result import HYPOTHESES, ColocResult


def _matrix(bf, snps, null_index, name):
    snps = val.ids(snps, name + " SNP IDs")
    a = np.asarray(bf, dtype=np.float64)
    if a.ndim == 1:
        a = a[None, :]
    n = len(snps)
    if null_index is not None and (isinstance(null_index, bool) or null_index != n):
        raise ValueError("an explicit null must be the last column, index len(snps)")
    if a.ndim != 2 or a.shape[0] == 0 or a.shape[1] != n + (null_index is not None):
        raise ValueError(f"{name} must have shape (L, m), or (L, m+1) with explicit null_index=m")
    if np.isnan(a).any() or np.isposinf(a).any():
        raise ValueError(f"{name} permits finite logBF or -inf, but not NaN or +inf")
    if null_index is not None:
        if not np.isfinite(a[:, -1]).all():
            raise ValueError("null logBF must be finite")
        a = a - a[:, -1, None]
    if not np.isfinite(a[:, :n]).any(axis=1).all():
        raise ValueError("every component must have some finite biological logBF")
    return a, snps


def _posterior_support(a, common_indices, n, p, has_null):
    # Conditioning on a biological association cancels the null mass exactly.
    # Compute this ratio on the full support, before SNP intersection.
    bio_norm = logsumexp(a[:, :n], axis=1)
    overlap = np.exp(logsumexp(a[:, common_indices], axis=1) - bio_norm)
    if not has_null:
        return overlap, np.zeros(len(a)), np.argmax(a[:, :n], axis=1)
    pi = min(p, 1 / n)
    prior = np.r_[np.repeat(pi, n), 1 - n * pi]
    prior[prior == 0] = 1e-16  # pinned logbf_to_pp zero-prior floor
    prior /= prior.sum()
    logpost = a + np.log(prior)
    norm = logsumexp(logpost, axis=1)
    return overlap, np.exp(logpost[:, -1] - norm), np.argmax(logpost, axis=1)


def coloc_bf_bf(
    bf1,
    bf2,
    *,
    snps1,
    snps2,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 5e-6,
    prior_weights1=None,
    prior_weights2=None,
    null_index1: int | None = None,
    null_index2: int | None = None,
    overlap_min: float = 0.5,
    trim_by_posterior: bool = True,
    engine: str = "numpy",
) -> ColocResult:
    """Colocalise every pair of natural-log Bayes-factor rows.

    Parameters
    ----------
    bf1, bf2 : array-like, shape (L, m) or (m,)
        Full component rows of natural log BF; -inf denotes zero BF. An
        explicit final null column changes shape to (L, m+1). No PIPs.
    snps1, snps2 : sequence of str
        Unique biological SNP IDs in each matrix's column order, excluding null.
    p1, p2, p12 : float
        Positive scalar priors, defaults 1e-4, 1e-4, 5e-6. BF hypothesis priors
        are not adjusted for region size by the upstream function.
    prior_weights1, prior_weights2 : mapping, Series or array, optional
        Nonnegative finite biological-SNP weights, bound by ID (arrays bind
        to supplied SNP order). Each side needs positive total on common SNPs.
        Unlike R's matrix-column length check, null weights are not requested.
    null_index1, null_index2 : int, optional
        Must equal m for an explicit final null. Null is never inferred.
    overlap_min : float
        Minimum fraction of each component's full biological posterior support
        represented by common SNPs. Coloc weights do not change this check.
    trim_by_posterior : bool
        Drop insufficient pairs. If all fail, retain reference-compatible rows
        with H1–H4 missing and H0=min(null posterior1, null posterior2).
    engine : {'numpy', 'numba'}
        O(m) stable H3 calculation without a SNP-by-SNP configuration matrix.

    Returns
    -------
    ColocResult
        First trait row varies fastest, matching R expand.grid. idx1/idx2 are
        zero based, idx1_r/idx2_r are R indices. SNP.PP.H4.abf (one pair) or
        SNP.PP.H4.rowK (multiple) gives P(shared SNP | H4, data), including
        both coloc weights. Parallel .weighted fields are identical aliases.

    Raises
    ------
    ValueError
        Invalid BF, identifiers, priors, weights or null convention.

    Notes
    -----
    Common columns and weights are always matched by SNP ID. No-overlap,
    insufficient-overlap and partial trimming are distinct result statuses.

    References
    ----------
    coloc 6.0.3 susie.R: coloc.bf_bf and logbf_to_pp.

    Examples
    --------
    >>> r = coloc_bf_bf([1., 3.], [2., 4.], snps1=['a', 'b'], snps2=['a', 'b'])
    >>> r.summary.idx1.tolist()
    [0]
    """
    pri = val.priors(p1, p2, p12)
    cutoff = val.scalar(overlap_min, "overlap_min", positive=False)
    if not 0 <= cutoff <= 1 or not isinstance(trim_by_posterior, (bool, np.bool_)):
        raise ValueError("overlap_min must be in [0,1] and trim_by_posterior must be bool")
    a, s1 = _matrix(bf1, snps1, null_index1, "bf1")
    b, s2 = _matrix(bf2, snps2, null_index2, "bf2")
    w1 = val.weights(prior_weights1, s1, "prior_weights1")
    w2 = val.weights(prior_weights2, s2, "prior_weights2")
    index2 = {s: j for j, s in enumerate(s2)}
    common = np.array([s for s in s1 if s in index2], dtype=str)
    meta = {
        "input_nsnps": [len(s1), len(s2)],
        "null_indices": [null_index1, null_index2],
        "overlap_min": cutoff,
        "trim_by_posterior": bool(trim_by_posterior),
    }
    if not len(common):
        return ColocResult(
            pd.DataFrame(columns=["nsnps", *HYPOTHESES, "idx1", "idx2"]),
            pd.DataFrame({"snp": common}),
            pri,
            "coloc.bf_bf",
            "no_overlap",
            metadata=meta,
        )
    index1 = {s: j for j, s in enumerate(s1)}
    i1, i2 = [index1[s] for s in common], [index2[s] for s in common]
    prop1, ph01, hit1 = _posterior_support(a, i1, len(s1), p1, null_index1 is not None)
    prop2, ph02, hit2 = _posterior_support(b, i2, len(s2), p2, null_index2 is not None)
    todo = [(i, j) for j in range(len(b)) for i in range(len(a))]
    keep = [
        (i, j)
        for i, j in todo
        if not trim_by_posterior or (prop1[i] >= cutoff and prop2[j] >= cutoff)
    ]
    meta["dropped_pairs"] = [list(pair) for pair in todo if pair not in keep]
    meta["overlap1"], meta["overlap2"] = prop1.tolist(), prop2.tolist()
    w1, w2 = val.subset_weights(w1, common), val.subset_weights(w2, common)
    weights = {"prior_weights1": w1, "prior_weights2": w2}
    if not keep:
        warnings.warn(
            "common SNPs carry insufficient posterior support",
            UserWarning,
            stacklevel=2,
        )
        names1, names2 = list(s1) + ["null"], list(s2) + ["null"]
        rows = [
            dict(
                nsnps=len(common),
                hit1=names1[hit1[i]],
                hit2=names2[hit2[j]],
                idx1=i,
                idx2=j,
                idx1_r=i + 1,
                idx2_r=j + 1,
                overlap1=prop1[i],
                overlap2=prop2[j],
                status="insufficient_overlap",
                **dict(zip(HYPOTHESES, [min(ph01[i], ph02[j]), *[np.nan] * 4])),
            )
            for i, j in todo
        ]
        return ColocResult(
            pd.DataFrame(rows),
            pd.DataFrame({"snp": common}),
            pri,
            "coloc.bf_bf",
            "insufficient_overlap",
            weights,
            meta,
        )
    a, b = a[:, i1], b[:, i2]
    rows, columns = [], {"snp": common}
    for k, (i, j) in enumerate(keep):
        pp, spp, log_u = combine(a[i], b[j], pri, w1, w2, engine)
        column = "SNP.PP.H4.abf" if len(keep) == 1 else f"SNP.PP.H4.row{k + 1}"
        columns[column], columns[column + ".weighted"] = spp, spp
        rows.append(
            dict(
                nsnps=len(common),
                hit1=common[np.argmax(a[i])],
                hit2=common[np.argmax(b[j])],
                idx1=i,
                idx2=j,
                idx1_r=i + 1,
                idx2_r=j + 1,
                overlap1=prop1[i],
                overlap2=prop2[j],
                status="ok",
                logU3=log_u[3],
                logU4=log_u[4],
                h4_status=h4_status(log_u[4], pp[4]),
                snp_pp_column=column,
                **dict(zip(HYPOTHESES, pp)),
            )
        )
    status = "ok" if len(keep) == len(todo) else "partially_trimmed"
    return ColocResult(
        pd.DataFrame(rows), pd.DataFrame(columns), pri, "coloc.bf_bf", status, weights, meta
    )


def _fit_dict(fit):
    if isinstance(fit, Mapping):
        return fit
    if hasattr(fit, "as_dict"):
        return fit.as_dict()
    if hasattr(fit, "__dict__"):
        return vars(fit)
    raise TypeError("fit must be a mapping or an object with as_dict()")


def _fit(fit):
    d = _fit_dict(fit)
    required = {"variant_ids", "lbf_variable", "sets", "converged"}
    if not required <= d.keys():
        raise ValueError(f"fit is missing {sorted(required - d.keys())}")
    a, snps = _matrix(d["lbf_variable"], d["variant_ids"], d.get("null_index"), "lbf_variable")
    sets = d["sets"]
    if not isinstance(sets, Mapping):
        raise ValueError("sets must be a mapping with cs and cs_index")
    cs = sets.get("cs", [])
    if cs is None:
        cs = []
    idx = np.asarray(sets.get("cs_index", []))
    if (
        idx.ndim != 1
        or len(idx) != len(cs)
        or (
            len(idx)
            and (
                not np.issubdtype(idx.dtype, np.integer)
                or len(set(idx)) != len(idx)
                or np.any(idx < 0)
                or np.any(idx >= len(a))
            )
        )
    ):
        raise ValueError(
            "cs_index must contain distinct valid zero-based component IDs, one per CS"
        )
    members = list(cs.values()) if isinstance(cs, Mapping) else list(cs)
    for member in members:
        x = np.asarray(member)
        if (
            x.ndim != 1
            or not len(x)
            or not np.issubdtype(x.dtype, np.integer)
            or np.any(x < 0)
            or np.any(x >= len(snps))
        ):
            raise ValueError("credible sets require nonempty biological SNP indices")
    if not isinstance(d["converged"], (bool, np.bool_)):
        raise ValueError("converged must be a boolean")
    return d, a, snps, idx.astype(np.int64), members


def coloc_susie(
    fit1,
    fit2,
    *,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 5e-6,
    prior_weights1=None,
    prior_weights2=None,
    overlap_min: float = 0.5,
    trim_by_posterior: bool = True,
    engine: str = "numpy",
) -> ColocResult:
    """Consume existing SuSiE fits and pair rows selected by actual cs_index.

    Parameters
    ----------
    fit1, fit2 : mapping or object with as_dict()
        Fields: variant_ids (m biological IDs), lbf_variable (L,m or L,m+1),
        sets.cs, sets.cs_index (zero-based original component IDs), converged,
        and optional null_index=m. Every SNP column is retained, including
        SNPs outside credible sets. Fine-mapping priors are not reapplied.
    p1, p2, p12 : float
        Positive scalar coloc priors, defaults 1e-4, 1e-4, 5e-6.
    prior_weights1, prior_weights2 : mapping, Series or array, optional
        Nonnegative finite coloc weights on biological SNPs of each fit,
        with positive total on each side after intersection.
    overlap_min : float
        Posterior overlap threshold, default 0.5.
    trim_by_posterior : bool
        Apply full-support overlap trimming, default True.
    engine : {'numpy', 'numba'}
        Stable hypothesis kernel implementation.

    Returns
    -------
    ColocResult
        Pair-specific hypotheses. Original zero-based component IDs and R
        index mapping are explicit. No CS returns status='no_cs'.

    Raises
    ------
    ValueError
        Invalid result contract or input, including invalid CS component IDs.

    Notes
    -----
    Unconverged fits produce a warning and metadata, preserving computable
    reference behavior. No fitting, R calls, or prusie import occurs here.

    References
    ----------
    coloc 6.0.3 susie.R: coloc.susie and coloc.bf_bf.

    Examples
    --------
    Pass the SusieResult returned by prusie.susie_rss directly as either fit.
    """
    pri = val.priors(p1, p2, p12)
    d1, a, s1, idx1, cs1 = _fit(fit1)
    d2, b, s2, idx2, cs2 = _fit(fit2)
    val.weights(prior_weights1, s1, "prior_weights1")
    val.weights(prior_weights2, s2, "prior_weights2")
    converged = [bool(d1["converged"]), bool(d2["converged"])]
    if not all(converged):
        warnings.warn(
            "at least one SuSiE fit did not converge; inspect fit diagnostics",
            UserWarning,
            stacklevel=2,
        )
    meta = {
        "converged": converged,
        "credible_sets1": {
            int(i): s1[np.asarray(c, dtype=int)].tolist() for i, c in zip(idx1, cs1)
        },
        "credible_sets2": {
            int(i): s2[np.asarray(c, dtype=int)].tolist() for i, c in zip(idx2, cs2)
        },
    }
    if not len(idx1) or not len(idx2):
        warnings.warn("at least one fit has no credible sets", UserWarning, stacklevel=2)
        return ColocResult(
            pd.DataFrame(columns=["nsnps", *HYPOTHESES, "idx1", "idx2"]),
            pd.DataFrame(columns=["snp"]),
            pri,
            "coloc.susie",
            "no_cs",
            metadata=meta,
        )
    result = coloc_bf_bf(
        a[idx1],
        b[idx2],
        snps1=s1,
        snps2=s2,
        **pri,
        prior_weights1=prior_weights1,
        prior_weights2=prior_weights2,
        null_index1=d1.get("null_index"),
        null_index2=d2.get("null_index"),
        overlap_min=overlap_min,
        trim_by_posterior=trim_by_posterior,
        engine=engine,
    )
    if len(result.summary):
        original1 = idx1[result.summary.idx1.to_numpy(dtype=int)]
        original2 = idx2[result.summary.idx2.to_numpy(dtype=int)]
        result.summary["idx1"] = original1
        result.summary["idx2"] = original2
        result.summary["idx1_r"] = original1 + 1
        result.summary["idx2_r"] = original2 + 1
    if "dropped_pairs" in result.metadata:
        result.metadata["dropped_pairs"] = [
            [int(idx1[i]), int(idx2[j])] for i, j in result.metadata["dropped_pairs"]
        ]
    result.method = "coloc.susie"
    result.metadata.update(meta)
    return result
