"""Simple tabular results; component indices are always zero based."""

from dataclasses import dataclass, field
from typing import Any

import pandas as pd

HYPOTHESES = [f"PP.H{i}.abf" for i in range(5)]
RESULT_SCHEMA_VERSION = "2"


@dataclass
class ColocResult:
    """Colocalisation output.

    Parameters
    ----------
    summary : pandas.DataFrame
        One row per comparison, with nsnps and PP.H0.abf through PP.H4.abf.
        Pairing adds idx1/idx2 (original zero-based component IDs),
        idx1_r/idx2_r (R IDs), hit1/hit2, overlap1/overlap2 and status.
        Computed rows include logU3, logU4 and h4_status: defined,
        zero_shared_evidence, or posterior_underflow. Only zero_shared_evidence
        makes the conditional posterior undefined; finite logU4 retains it.
    results : pandas.DataFrame
        SNP IDs and P(shared SNP | H4, data), using both coloc weights in the
        main SNP.PP.H4 fields. Parallel .weighted fields are equal aliases.
        These are conditional probabilities, not P(H4 and shared SNP | data).
    priors : dict
        Actual scalar priors after any ABF adjustment.
    method : str
        coloc.abf, coloc.bf_bf or coloc.susie.
    status : str
        ok, partially_trimmed, insufficient_overlap, no_overlap, or no_cs.
    weights : dict
        Common-SNP indexed weight series, or None for an unweighted side.
    metadata : dict
        Diagnostics, original support, convergence and selected credible sets.

    Notes
    -----
    There is no aggregate H4 for multiple signal pairs. Empty results do not
    imply H4=0. ``result['summary']`` and attribute access are equivalent.
    """

    summary: pd.DataFrame
    results: pd.DataFrame
    priors: dict[str, float]
    method: str
    status: str = "ok"
    weights: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.metadata.setdefault("result_schema_version", RESULT_SCHEMA_VERSION)

    @property
    def nsnps(self) -> int:
        """Number of common biological SNPs, excluding an explicit null."""
        return int(self.summary.iloc[0]["nsnps"]) if len(self.summary) else 0

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def as_dict(self) -> dict[str, Any]:
        """Return a shallow view without requiring a third-party result class."""
        return vars(self).copy()
