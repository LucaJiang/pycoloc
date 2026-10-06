"""Wakefield approximate Bayes factors and single-causal-variant coloc."""

import warnings
from collections.abc import Mapping

import numpy as np
import pandas as pd
from scipy.special import ndtri_exp

from . import _validation as val
from ._math import combine, h4_status
from .result import HYPOTHESES, ColocResult


def process_dataset(dataset: Mapping, effect_priors: Mapping | None = None):
    if not isinstance(dataset, Mapping):
        raise TypeError("dataset must be a mapping")
    if dataset.get("type") not in ("quant", "cc"):
        raise ValueError("type must be 'quant' or 'cc'")
    if "snp" not in dataset:
        raise ValueError("dataset requires explicit snp IDs")
    snps = val.ids(dataset["snp"])
    n = len(snps)
    arrays = {
        k: val.vector(dataset[k], n, k)
        for k in ("beta", "varbeta", "pvalues", "MAF", "position")
        if k in dataset
    }
    if "varbeta" in arrays and np.any(arrays["varbeta"] <= 0):
        raise ValueError("varbeta must be positive; square SE before calling")
    if "MAF" in arrays and np.any((arrays["MAF"] <= 0) | (arrays["MAF"] >= 1)):
        raise ValueError("MAF must lie strictly between zero and one")
    if "pvalues" in arrays and np.any((arrays["pvalues"] < 0) | (arrays["pvalues"] > 1)):
        raise ValueError("pvalues must lie between zero and one")
    N = val.scalar(dataset["N"], "N") if "N" in dataset else None
    s = val.scalar(dataset["s"], "s") if "s" in dataset else None
    if s is not None and s >= 1:
        raise ValueError("case fraction s must lie strictly between zero and one")
    sd_y = val.scalar(dataset["sdY"], "sdY") if "sdY" in dataset else None
    prior = {"quant": 0.15, "cc": 0.2}
    if effect_priors is not None:
        if not isinstance(effect_priors, Mapping) or set(effect_priors) != set(prior):
            raise ValueError("effect_priors must specify quant and cc standard deviations")
        prior = {k: val.scalar(v, f"effect_priors[{k}]") for k, v in effect_priors.items()}
    estimated = False
    if "beta" in arrays and "varbeta" in arrays:
        V = arrays["varbeta"]
        z = arrays["beta"] / np.sqrt(V)
        if dataset["type"] == "quant" and sd_y is None:
            if N is None or "MAF" not in arrays:
                raise ValueError("quant beta/varbeta requires sdY or MAF and N")
            warnings.warn(
                "estimating sdY from MAF and varbeta; supply known sdY when available",
                UserWarning,
                stacklevel=3,
            )
            # No-intercept least squares: nvx ~ 1/varbeta. Scale the design
            # before taking inner products to avoid squaring huge precisions.
            vmin = np.min(V)
            x = vmin / V
            y = 2 * N * arrays["MAF"] * (1 - arrays["MAF"])
            sd_y = np.sqrt(vmin * np.dot(x, y) / np.dot(x, x))
            estimated = True
        sd_prior = prior["quant"] * sd_y if dataset["type"] == "quant" else prior["cc"]
        branch = "beta/varbeta"
    else:
        if effect_priors is not None:
            raise ValueError(
                "effect_priors is unsupported for pvalues/MAF/N; "
                "omit it for fixed prior SDs or supply beta/varbeta"
            )
        if not all(k in arrays for k in ("pvalues", "MAF")) or N is None:
            raise ValueError("supply beta/varbeta or pvalues/MAF/N")
        if np.any(arrays["pvalues"] == 0):
            raise ValueError(
                "pvalue ABF requires positive P; supply beta/varbeta for underflowed P"
            )
        f = arrays["MAF"]
        V = 1 / (2 * N * f * (1 - f))
        if dataset["type"] == "cc":
            if s is None:
                raise ValueError("case-control pvalue ABF requires case fraction s")
            V = V / (s * (1 - s))
        z = -ndtri_exp(np.log(arrays["pvalues"]) - np.log(2.0))
        sd_prior = 0.15 if dataset["type"] == "quant" else 0.2
        branch = "pvalues/MAF/N"
    log_ratio = 2 * np.log(sd_prior) - np.log(V)
    log_shrink = -np.logaddexp(0, -log_ratio)
    r = np.exp(log_shrink)
    labf = 0.5 * (-np.logaddexp(0, log_ratio) + r * z**2)
    if not np.isfinite(labf).all():
        raise ValueError("ABFs exceed float64 range; check effect and variance units")
    df = pd.DataFrame({"snp": snps, "V": V, "z": z, "r": r, "lABF": labf})
    for k in ("position", "pvalues", "MAF"):
        if k in arrays:
            df[k] = arrays[k]
    return df, {
        "branch": branch,
        "sdY": sd_y,
        "sdY_estimated": estimated,
        "effect_prior_sd": sd_prior,
    }


def coloc_abf(
    dataset1: Mapping,
    dataset2: Mapping,
    *,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
    prior_weights1=None,
    prior_weights2=None,
    MAF=None,
    effect_priors: Mapping | None = None,
    engine: str = "numpy",
) -> ColocResult:
    """Colocalise two traits under at most one causal variant per trait.

    Parameters
    ----------
    dataset1, dataset2 : mapping
        Unique ``snp`` IDs, ``type`` ('quant' or 'cc'), and either beta/varbeta
        or pvalues/MAF/N vectors. beta for cc is log(OR); varbeta is SE².
        Quant beta requires sdY or MAF/N for official sdY estimation. The
        pvalue cc path also requires scalar case fraction s. beta/varbeta
        takes precedence. N is a scalar; no LD is required.
    p1, p2, p12 : float
        Positive scalar SNP priors. Defaults 1e-4, 1e-4, 1e-5; large-region
        adjustment follows pinned coloc 6.0.3 before/after intersection.
    prior_weights1, prior_weights2 : mapping, pandas.Series or array, optional
        Nonnegative finite weights, positive total on common SNPs. Keys must
        exactly match each input's SNP IDs; arrays bind to input SNP order.
        Intersection is followed by renormalisation. Fine-mapping priors are
        a separate parameter and must not be copied automatically here.
    MAF : mapping or array, optional
        Shared fallback frequency. Arrays require identical input SNP order.
    effect_priors : mapping, optional
        Prior SDs, default {'quant': 0.15, 'cc': 0.2}. Quant is multiplied by
        sdY on the estimate path. Explicit nondefaults are honored here,
        unlike the pinned R public wrapper which silently drops this option.
        P-value ABFs use fixed SDs (quant: 0.15 on a unit-variance phenotype;
        cc: 0.2 on the logOR scale). Explicit effect_priors with an active
        P-value branch raises ValueError. SDs are squared to obtain variance W.
    engine : {'numpy', 'numba'}
        Stable O(m) H3 reduction. Numba is optional, without fastmath.

    Returns
    -------
    ColocResult
        One-row summary, per-SNP ABFs and conditional shared probabilities.

    Raises
    ------
    ValueError
        Invalid statistics, priors or weights. No overlap instead returns
        status='no_overlap' with an empty summary, never H4=0.

    Notes
    -----
    SNP.PP.H4 is conditional on H4 and includes both coloc weights;
    SNP.PP.H4.weighted is an alias. SNPs are aligned by ID; conflicting
    shared positions raise rather than silently changing the analysis support.

    References
    ----------
    coloc 6.0.3, claudia.R: process.dataset, combine_abf_weighted, coloc.abf.

    Examples
    --------
    >>> d = dict(snp=['a', 'b'], beta=[0.2, 0.0], varbeta=[0.01, 0.01], type='cc')
    >>> coloc_abf(d, d).nsnps
    2
    """
    pri = val.priors(p1, p2, p12)
    d1, d2 = dict(dataset1), dict(dataset2)
    if MAF is not None:
        if not isinstance(MAF, Mapping) and not np.array_equal(d1["snp"], d2["snp"]):
            raise ValueError("shared MAF array requires identical SNP order; use keyed MAF")
        for d in (d1, d2):
            if "MAF" not in d:
                d["MAF"] = [MAF[s] for s in d["snp"]] if isinstance(MAF, Mapping) else MAF
    df1, meta1 = process_dataset(d1, effect_priors)
    df2, meta2 = process_dataset(d2, effect_priors)
    w1 = val.weights(prior_weights1, df1.snp.to_numpy(), "prior_weights1")
    w2 = val.weights(prior_weights2, df2.snp.to_numpy(), "prior_weights2")
    if w1 is not None and w2 is not None:
        warnings.warn(
            "both trait weights supplied; derive them from independent information",
            UserWarning,
            stacklevel=2,
        )
    pri["p1"] = val.adjust_prior(pri["p1"], len(df1), "p1")
    pri["p2"] = val.adjust_prior(pri["p2"], len(df2), "p2")
    merged = df1.merge(df2, on="snp", suffixes=(".df1", ".df2"), sort=True)
    meta = {"dataset1": meta1, "dataset2": meta2, "input_nsnps": [len(df1), len(df2)]}
    if not len(merged):
        return ColocResult(
            pd.DataFrame(columns=["nsnps", *HYPOTHESES]),
            merged,
            pri,
            "coloc.abf",
            "no_overlap",
            metadata=meta,
        )
    if "position.df1" in merged and "position.df2" in merged:
        if not np.array_equal(merged["position.df1"], merged["position.df2"]):
            raise ValueError("shared SNPs have conflicting positions")
        merged["position"] = merged.pop("position.df1")
        merged.pop("position.df2")
    pri["p12"] = val.adjust_prior(pri["p12"], len(merged), "p12")
    w1, w2 = val.subset_weights(w1, merged.snp), val.subset_weights(w2, merged.snp)
    l1, l2 = merged["lABF.df1"].to_numpy(), merged["lABF.df2"].to_numpy()
    pp, snp_pp, log_u = combine(l1, l2, pri, w1, w2, engine)
    merged["internal.sum.lABF"] = l1 + l2
    merged["SNP.PP.H4"], merged["SNP.PP.H4.weighted"] = snp_pp, snp_pp
    summary = pd.DataFrame(
        [
            dict(
                nsnps=len(merged),
                logU3=log_u[3],
                logU4=log_u[4],
                h4_status=h4_status(log_u[4], pp[4]),
                **dict(zip(HYPOTHESES, pp)),
            )
        ]
    )
    return ColocResult(
        summary,
        merged,
        pri,
        "coloc.abf",
        weights={"prior_weights1": w1, "prior_weights2": w2},
        metadata=meta,
    )
