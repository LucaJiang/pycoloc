"""Identity-bound LD data for plotting; no LD estimation or allele recoding."""

from dataclasses import dataclass, field

import numpy as np


def _ids(value):
    ids = np.asarray(value, dtype=str).copy()
    if ids.ndim != 1 or not len(ids) or len(set(ids)) != len(ids) or np.any(ids == ""):
        raise ValueError("LD requires nonempty unique SNP IDs")
    return ids


def _r2(value, kind, tolerance):
    if kind not in {"r", "r2"}:
        raise ValueError("LD kind must explicitly be 'r' or 'r2'")
    if not np.isfinite(tolerance) or not 0 <= tolerance <= 1e-4:
        raise ValueError("LD tolerance must lie in [0, 1e-4]")
    a = np.asarray(value, dtype=float)
    low = -1 if kind == "r" else 0
    if not np.isfinite(a).all() or np.any(a < low - tolerance) or np.any(a > 1 + tolerance):
        raise ValueError("LD must be finite and within its declared range")
    a = np.clip(a, low, 1)
    return a * a if kind == "r" else a.copy()


@dataclass
class LDMatrix:
    """Pairwise LD with explicit ID/coordinate order and signed-r or r² units.

    Symmetry/range/diagonal tolerance defaults to 1e-8. Only tolerated roundoff
    is clipped; matrices are never silently symmetrized. Inputs are copied.
    Provenance should identify population, samples, estimator and source files.
    """

    values: object
    snp_ids: object
    positions: object
    kind: str
    panel: str
    provenance: dict = field(default_factory=dict)
    tolerance: float = 1e-8

    def __post_init__(self):
        self.snp_ids = _ids(self.snp_ids)
        self.values = np.asarray(self.values, dtype=float).copy()
        n = len(self.snp_ids)
        self.positions = np.asarray(self.positions, dtype=float).copy()
        if self.values.shape != (n, n):
            raise ValueError("LD matrix shape does not match SNP IDs")
        if (
            self.positions.shape != (n,)
            or not np.isfinite(self.positions).all()
            or np.any(self.positions < 1)
            or np.any(self.positions != np.floor(self.positions))
        ):
            raise ValueError("LD positions must be positive integer base-pair coordinates")
        self.r2 = _r2(self.values, self.kind, self.tolerance)
        if not np.allclose(self.values, self.values.T, atol=self.tolerance, rtol=0):
            raise ValueError("LD matrix must be symmetric")
        if not np.allclose(np.diag(self.values), 1, atol=self.tolerance, rtol=0):
            raise ValueError("LD diagonal must equal one")
        if not self.panel:
            raise ValueError("An explicit LD panel label is required")


@dataclass
class LeadLD:
    """An ID-bound lead vector for association colors, never pairwise heatmaps."""

    values: object
    snp_ids: object
    lead_id: str
    kind: str
    panel: str
    provenance: dict = field(default_factory=dict)
    tolerance: float = 1e-8

    def __post_init__(self):
        self.snp_ids = _ids(self.snp_ids)
        self.r2 = _r2(self.values, self.kind, self.tolerance)
        if self.r2.shape != self.snp_ids.shape or not self.lead_id or not self.panel:
            raise ValueError("Lead LD requires matching vector/IDs, reference lead and panel")
        index = np.flatnonzero(self.snp_ids == self.lead_id)
        if len(index) and abs(self.r2[index[0]] - 1) > self.tolerance * 2:
            raise ValueError("Reference lead LD with itself must equal one")


def as_ld(value):
    if value is None or isinstance(value, (LDMatrix, LeadLD)):
        return value
    if isinstance(value, dict):
        return (LeadLD if "lead_id" in value else LDMatrix)(**value)
    raise TypeError("LD must be LDMatrix, LeadLD, an equivalent mapping, or None")


def lead_colors(ld, ids, positions, lead):
    """Map r² by identity; missing or unknown reference leads stay NaN."""
    result = np.full(len(ids), np.nan)
    if ld is None:
        return result, "LD not supplied"
    mapping = dict(zip(ld.snp_ids, range(len(ld.snp_ids))))
    if isinstance(ld, LDMatrix):
        for snp, position in zip(ids, positions):
            if snp in mapping and ld.positions[mapping[snp]] != position:
                raise ValueError("LD and association coordinates disagree for " + snp)
        if lead not in mapping:
            return result, "Selected lead absent from LD panel"
        values = ld.r2[mapping[lead]]
    else:
        if lead != ld.lead_id:
            raise ValueError("Lead vector reference differs from selected association lead")
        values = ld.r2
    for j, snp in enumerate(ids):
        if snp in mapping:
            result[j] = values[mapping[snp]]
    return result, "ID-aligned lead r²"


def heatmap_selection(ld, start, end, max_snps, lead, allowed_ids=None):
    if not isinstance(ld, LDMatrix):
        raise ValueError("A pairwise LDMatrix is required for a heatmap")
    if not isinstance(max_snps, int) or max_snps < 2:
        raise ValueError("max_snps must be an integer of at least two")
    order = np.lexsort((ld.snp_ids, ld.positions))
    order = order[(ld.positions[order] >= start) & (ld.positions[order] <= end)]
    allowed = None if allowed_ids is None else set(allowed_ids)
    if allowed is not None:
        order = np.asarray([j for j in order if ld.snp_ids[j] in allowed], dtype=int)
    # Repeated physical positions cannot define distinct physical cells. Keep a
    # deterministic representative, preferring the requested lead, and audit it.
    unique = []
    omitted_duplicates = []
    for position in np.unique(ld.positions[order]):
        group = order[ld.positions[order] == position]
        chosen = next((int(j) for j in group if ld.snp_ids[j] == lead), int(group[0]))
        unique.append(chosen)
        omitted_duplicates.extend(int(j) for j in group if j != chosen)
    indices = np.asarray(unique, dtype=int)
    if len(indices) > max_snps:
        keep = np.linspace(0, len(indices) - 1, max_snps, dtype=int)
        lead_rank = np.flatnonzero(ld.snp_ids[indices] == lead)
        if len(lead_rank) and lead_rank[0] not in keep:
            replace = int(np.argmin(np.abs(keep - lead_rank[0])))
            keep[replace] = lead_rank[0]
            keep.sort()
        indices = indices[keep]
    selected = set(indices.tolist())
    metadata = dict(
        panel=ld.panel,
        provenance=ld.provenance,
        input_kind=ld.kind,
        tolerance=ld.tolerance,
        geometry="physical bp; x=(u+v)/2, y=(v-u)/2 for cell edges",
        rule="sort position then SNP ID; one representative per position (lead preferred); evenly spaced ranks; preserve lead",
        selected_ids=ld.snp_ids[indices].tolist(),
        positions=ld.positions[indices].astype(int).tolist(),
        original_indices=indices.tolist(),
        input_count=len(ld.snp_ids),
        max_snps=max_snps,
        omitted_ids=[str(ld.snp_ids[j]) for j in range(len(ld.snp_ids)) if j not in selected],
        repeated_position_omissions=ld.snp_ids[omitted_duplicates].tolist(),
        lead=lead,
        lead_retained=lead in set(ld.snp_ids[indices]),
        omitted_outside_region=ld.snp_ids[(ld.positions < start) | (ld.positions > end)].tolist(),
        allowed_id_count=None if allowed is None else len(allowed),
        omitted_not_shared=[]
        if allowed is None
        else [str(s) for s in ld.snp_ids if s not in allowed],
    )
    return indices, metadata
