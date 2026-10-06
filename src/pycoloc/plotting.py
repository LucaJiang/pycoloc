"""Editable regional marginal-association tracks with result-derived labels."""

from collections.abc import Mapping

import numpy as np
from scipy.special import log_ndtr

from . import _validation as val
from .result import HYPOTHESES, ColocResult


def _track(dataset: Mapping):
    if "snp" not in dataset or "position" not in dataset or "chromosome" not in dataset:
        raise ValueError("plot datasets require snp, position and chromosome")
    snps = val.ids(dataset["snp"])
    n = len(snps)
    pos = val.vector(dataset["position"], n, "position")
    if np.any(pos <= 0):
        raise ValueError("positions must be positive base-pair coordinates")
    chrom = np.asarray(dataset["chromosome"])
    if chrom.ndim == 0:
        chrom = np.repeat(str(chrom), n)
    if chrom.shape != (n,) or len(set(chrom.astype(str))) != 1:
        raise ValueError("one chromosome is required per regional plot")
    source = ""
    z = None
    if "beta" in dataset and "varbeta" in dataset:
        v = val.vector(dataset["varbeta"], n, "varbeta")
        if np.any(v <= 0):
            raise ValueError("varbeta must be positive")
        z = val.vector(dataset["beta"], n, "beta") / np.sqrt(v)
    elif "z" in dataset:
        z = val.vector(dataset["z"], n, "z")
    logp = None
    if "logp" in dataset:
        logp = val.vector(dataset["logp"], n, "logp")
        if np.any(logp > 0):
            raise ValueError("logp must be natural logarithm of P, at most zero")
    if "pvalues" in dataset:
        p = val.vector(dataset["pvalues"], n, "pvalues")
        if np.any((p < 0) | (p > 1)):
            raise ValueError("plot P values must lie in [0,1]")
        positive = p > 0
        y = np.empty(n)
        y[positive] = -np.log10(p[positive])
        source = "original summary P"
        zero = ~positive
        if zero.any():
            if logp is not None:
                y[zero] = -logp[zero] / np.log(10)
                source += "; P=0 recovered from logP"
            elif z is not None:
                y[zero] = -(np.log(2) + log_ndtr(-np.abs(z[zero]))) / np.log(10)
                source += "; P=0 recovered from beta/SE or Z"
            else:
                y[zero] = -np.log10(np.nextafter(0.0, 1.0))
                source += "; P=0 clipped at smallest float64 (unrecoverable)"
    elif logp is not None:
        y = -logp / np.log(10)
        source = "supplied natural logP"
    elif z is not None:
        y = -(np.log(2) + log_ndtr(-np.abs(z))) / np.log(10)
        source = "two-sided Wald P from beta/SE or Z"
    else:
        raise ValueError("plot requires original pvalues, natural logp, or beta/varbeta or z")
    if not np.isfinite(y).all():
        raise ValueError("plot -log10(P) exceeds float64 range")
    return snps, pos, str(chrom[0]), y, source


def plot_coloc(
    result: ColocResult,
    dataset1: Mapping,
    dataset2: Mapping,
    *,
    signal_pair: tuple[int, int] | None = None,
    trait1: str = "Trait 1",
    trait2: str = "Trait 2",
    figsize: tuple[float, float] = (10, 7),
    fontsize: float = 10,
    point_size: float = 14,
    annotate_leads: bool = True,
    highlight_cs: bool = True,
):
    """Plot two regional marginal-association tracks and exact result summaries.

    Parameters
    ----------
    result : ColocResult
        Output of any coloc interface. Empty/ineligible states show no invented
        probabilities. Multiple valid pairs require explicit signal_pair.
    dataset1, dataset2 : mapping
        snp, position (bp), chromosome (scalar or constant vector), and pvalues,
        natural logp, beta/varbeta or z. Valid original P values take precedence.
        P=0 is recovered from logp or Wald statistics if available, otherwise
        explicitly clipped for display. Inputs are never changed.
    signal_pair : tuple of int, optional
        Original zero-based component IDs, not credible-set ordinals.
    trait1, trait2 : str
        Trait/gene names printed above the tracks.
    figsize : tuple of float
        Figure dimensions in inches.
    fontsize, point_size : float
        Local font size and marker area; no global matplotlib settings change.
    annotate_leads, highlight_cs : bool
        Mark actual result lead variants, the shared candidate selected from
        the main conditional SNP posterior, and selected component CS members.

    Returns
    -------
    figure, axes
        Editable matplotlib figure and two axes; never calls show().

    Raises
    ------
    ValueError
        Missing coordinates/statistics, chromosome conflicts, unknown signal
        pair, or plotted datasets missing analyzed SNPs.

    Notes
    -----
    These are marginal tracks, not conditional or residual associations. H0–H4
    come directly from a single result row and are never aggregated over pairs.
    Save with figure.savefig('region.svg') or figure.savefig('region.png').

    References
    ----------
    The statistical rows use the coloc 6.0.3 conventions documented by pycoloc.

    Examples
    --------
    >>> # fig, axes = plot_coloc(result, d1, d2, signal_pair=(0, 2))
    """
    import matplotlib.pyplot as plt

    tracks = [_track(dataset1), _track(dataset2)]
    if tracks[0][2] != tracks[1][2]:
        raise ValueError("both plot datasets must use the same chromosome")
    coordinates = [dict(zip(t[0], t[1])) for t in tracks]
    for snp in set(coordinates[0]) & set(coordinates[1]):
        if coordinates[0][snp] != coordinates[1][snp]:
            raise ValueError("shared SNPs have conflicting plot positions")
    analyzed = set(result.results.get("snp", []))
    for track in tracks:
        if not analyzed <= set(track[0]):
            raise ValueError("plot dataset is missing analyzed SNP IDs")
    if "position" in result.results:
        for snp, position in zip(result.results.snp, result.results.position):
            if coordinates[0][snp] != position:
                raise ValueError("plot positions disagree with analyzed positions")
    row = None
    summary = result.summary
    if signal_pair is not None:
        if "idx1" not in summary or "idx2" not in summary:
            raise ValueError("signal_pair only applies to BF/SuSiE pair results")
        selected = summary[(summary.idx1 == signal_pair[0]) & (summary.idx2 == signal_pair[1])]
        if len(selected) != 1:
            raise ValueError("signal_pair does not identify exactly one retained comparison")
        row = selected.iloc[0]
    elif len(summary) == 1:
        row = summary.iloc[0]
    elif len(summary) > 1:
        raise ValueError("multiple signal pairs: specify original component IDs with signal_pair")
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=figsize)
    fig.subplots_adjust(top=0.71, bottom=0.11, left=0.10, right=0.97, hspace=0.35)
    coords = np.concatenate([t[1] for t in tracks])
    region = f"chr{tracks[0][2]}:{coords.min():.0f}–{coords.max():.0f}"
    nsnps = int(row["nsnps"]) if row is not None else result.nsnps
    pair_text = ""
    if row is not None and "idx1" in row:
        pair_text = f" | components ({int(row.idx1)}, {int(row.idx2)}) [0-based]"
    title = f"{trait1} / {trait2} | {region}\n{result.method}{pair_text} | {nsnps} analysed SNPs"
    fig.text(0.5, 0.965, title, ha="center", va="top", fontsize=fontsize + 1)
    eligible = row is not None and all(np.isfinite(row[h]) for h in HYPOTHESES)
    probabilities = (
        "   ".join(f"PP.H{i} = {row[h]:.6g}" for i, h in enumerate(HYPOTHESES))
        if eligible
        else f"Status: {result.status}; no comparable H0–H4 posterior"
    )
    fig.text(0.5, 0.868, probabilities, ha="center", va="top", fontsize=fontsize)
    weighted = any(w is not None for w in result.weights.values())
    prior_text = ", ".join(f"{k}={v:.6g}" for k, v in result.priors.items())
    converged = result.metadata.get("converged", [True, True])
    if not all(converged):
        prior_text += " | WARNING: fit not converged"
    fig.text(
        0.5,
        0.817,
        prior_text + f" | variant-specific priors: {'yes' if weighted else 'no'}",
        ha="center",
        va="top",
        fontsize=fontsize - 1,
    )
    shared_candidate = None
    shared_column = None if row is None else row.get("snp_pp_column", "SNP.PP.H4")
    if eligible and shared_column in result.results:
        posterior = result.results[shared_column]
        if np.isfinite(posterior).any():
            shared_candidate = str(result.results.loc[posterior.idxmax(), "snp"])
    lead_ids = []
    for k, (ax, track, name) in enumerate(zip(axes, tracks, (trait1, trait2)), start=1):
        snps, pos, chrom, y, source = track
        x = pos / 1e6
        ax.scatter(
            x,
            y,
            s=point_size,
            color=("#3568a8" if k == 1 else "#b55c32"),
            linewidths=0,
            label="Marginal association",
        )
        lead = row.get(f"hit{k}") if row is not None else None
        if lead is None and result.method == "coloc.abf" and len(result.results):
            col = f"lABF.df{k}"
            lead = result.results.loc[result.results[col].idxmax(), "snp"]
        if highlight_cs and row is not None and f"idx{k}" in row:
            members = result.metadata.get(f"credible_sets{k}", {}).get(int(row[f"idx{k}"]), [])
            mask = np.isin(snps, members)
            if mask.any():
                ax.scatter(
                    x[mask],
                    y[mask],
                    s=point_size * 2,
                    facecolors="none",
                    edgecolors="#222222",
                    label="Selected CS",
                )
        if lead in set(snps):
            j = int(np.flatnonzero(snps == lead)[0])
            ax.scatter(
                [x[j]],
                [y[j]],
                marker="*",
                s=point_size * 7,
                c="#8b1b56",
                label="Component lead",
            )
            if annotate_leads:
                ax.annotate(
                    str(lead),
                    (x[j], y[j]),
                    xytext=(4, 6),
                    textcoords="offset points",
                    fontsize=fontsize - 2,
                )
            lead_ids.append(str(lead))
        if shared_candidate in set(snps):
            j = int(np.flatnonzero(snps == shared_candidate)[0])
            ax.scatter(
                [x[j]],
                [y[j]],
                marker="D",
                s=point_size * 3,
                facecolors="none",
                edgecolors="#26734d",
                label="Shared candidate | H4",
            )
        ax.set_title(
            f"{name} — marginal association\n{source}",
            loc="left",
            fontsize=fontsize - 1,
        )
        ax.set_ylabel("−log₁₀(P)", fontsize=fontsize)
        ax.tick_params(labelsize=fontsize - 1)
        ax.margins(x=0.03, y=0.20)
        ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel(f"Chromosome {tracks[0][2]} position (Mb)", fontsize=fontsize)
    fig.coloc_metadata = {
        "method": result.method,
        "region": region,
        "probabilities": probabilities,
        "signal_pair": None if row is None or "idx1" not in row else (int(row.idx1), int(row.idx2)),
        "leads": lead_ids,
        "shared_candidate": shared_candidate,
        "shared_posterior_column": shared_column,
        "h4_status": None if row is None else row.get("h4_status"),
        "p_sources": [t[4] for t in tracks],
    }
    return fig, axes
