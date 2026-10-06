"""Input and SNP-key validation shared by statistical interfaces."""

import warnings
from collections.abc import Mapping

import numpy as np
import pandas as pd


def ids(values, name="snp") -> np.ndarray:
    raw = np.asarray(values, dtype=object)
    if raw.ndim != 1 or not len(raw) or pd.isna(raw).any():
        raise ValueError(f"{name} must be a nonempty one-dimensional sequence without missing IDs")
    out = raw.astype(str)
    if np.any(out == "") or len(set(out)) != len(out):
        raise ValueError(f"{name} must contain unique nonempty IDs")
    return out


def vector(value, n: int, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=np.float64)
    if out.shape != (n,) or not np.isfinite(out).all():
        raise ValueError(f"{name} must be a finite vector of length {n}")
    return out


def scalar(value, name: str, *, positive=True) -> float:
    a = np.asarray(value, dtype=np.float64)
    if a.ndim != 0 or not np.isfinite(a) or (positive and a <= 0):
        raise ValueError(f"{name} must be a finite {'positive ' if positive else ''}scalar")
    return float(a)


def priors(p1, p2, p12) -> dict[str, float]:
    out = {k: scalar(v, k) for k, v in zip(("p1", "p2", "p12"), (p1, p2, p12))}
    if max(out.values()) > 1:
        raise ValueError("base priors must lie in (0, 1]; use weights for variant-specific priors")
    return out


def adjust_prior(p: float, n: int, name: str) -> float:
    if n * p >= 1:
        warnings.warn(
            f"{name} * nsnps >= 1; setting {name}=1/(nsnps+1)",
            UserWarning,
            stacklevel=3,
        )
        return 1.0 / (n + 1)
    return p


def weights(value, snps: np.ndarray, name: str) -> pd.Series | None:
    if value is None:
        return None
    if isinstance(value, (Mapping, pd.Series)):
        value = pd.Series(value, dtype=np.float64)
        value.index = value.index.astype(str)
        if not value.index.is_unique or set(value.index) != set(snps):
            raise ValueError(f"{name} keys must exactly match the biological SNP IDs")
        a = vector(value.reindex(snps).to_numpy(), len(snps), name)
    else:
        a = vector(value, len(snps), name)
    if np.any(a < 0) or not np.any(a > 0):
        raise ValueError(f"{name} must be nonnegative with positive total")
    return pd.Series(a, index=snps, name=name)


def subset_weights(w, snps):
    if w is None:
        return None
    out = w.reindex(snps)
    if not np.any(out.to_numpy() > 0):
        raise ValueError("weights have zero total on the common SNP support")
    return out
