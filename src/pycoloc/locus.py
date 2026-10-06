"""Five-panel, editable locus figures with offline annotation and explicit LD."""

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .annotation.models import (
    gene_exons,
    matches_gene,
    normalize_chromosome,
    select_transcript,
    union_exons,
    validate_region,
)
from .locus_ld import LDMatrix, as_ld, heatmap_selection, lead_colors
from .plotting import _track
from .result import HYPOTHESES

LD_BINS = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0])
LD_COLORS = ["#386CB0", "#56B4E9", "#009E73", "#E69F00", "#D55E00"]
MISSING_COLOR = "#B8B8B8"
LEAD_COLOR = "#762A83"


def _scale():
    from matplotlib.colors import BoundaryNorm, ListedColormap

    return ListedColormap(LD_COLORS), BoundaryNorm(LD_BINS, len(LD_COLORS), clip=True)


def _clean(value):
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return _clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _bound_text(text, ax, renderer, margin=3):
    """Translate text inside an axes rectangle; preserve the full label."""
    box = text.get_window_extent(renderer)
    bounds = ax.get_window_extent(renderer)
    if box.width > bounds.width - 2 * margin or box.height > bounds.height - 2 * margin:
        return False
    dx = max(bounds.x0 + margin - box.x0, 0) + min(bounds.x1 - margin - box.x1, 0)
    dy = max(bounds.y0 + margin - box.y0, 0) + min(bounds.y1 - margin - box.y1, 0)
    if dx or dy:
        transform = text.get_transform()
        old = transform.transform(text.get_position())
        text.set_position(transform.inverted().transform(old + np.array([dx, dy])))
    return True


def plot_gene_track(
    annotation, *, target_gene=None, ax=None, target_mode="union", max_lanes=6, fontsize=7
):
    """Draw all overlapping genes; default exons are the transcript union.

    ``target_mode='canonical'`` prefers MANE_Select, then Ensembl_canonical,
    then a clearly recorded representative fallback. Other genes remain unions.
    Labels are bounded and collision-filtered with target priority; gene bodies
    are retained even when their labels are suppressed. Full coordinates remain
    in ``ax.gene_metadata``; only visible segments are clipped to the window.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if target_mode not in {"union", "canonical"} or not isinstance(max_lanes, int) or max_lanes < 1:
        raise ValueError("target_mode is union/canonical and max_lanes must be positive")
    if ax is None:
        _, ax = plt.subplots(figsize=(7.09, 1.6))
    start, end = annotation.start, annotation.end
    targets = [g for g in annotation.genes if matches_gene(g, target_gene)]
    if target_gene is not None and len(targets) != 1:
        raise ValueError("Target gene must identify exactly one overlapping annotation gene")
    genes = sorted(
        annotation.genes, key=lambda g: (not matches_gene(g, target_gene), g.start, g.gene_id)
    )
    occupied = [[] for _ in range(max_lanes)]
    records, labels = [], []
    padding = (end - start + 1) * 0.016
    for gene in genes:
        target = matches_gene(gene, target_gene)
        left, right = max(start, gene.start), min(end, gene.end)
        scores = [
            sum(max(0, min(right + padding, b) - max(left - padding, a)) for a, b in intervals)
            for intervals in occupied
        ]
        lane = int(np.argmin(scores))
        occupied[lane].append((left - padding, right + padding))
        color = "#222222" if target else "#92999F"
        transcript, rule = (None, "collapsed union of all transcript exons")
        exons = gene_exons(gene)
        if target and target_mode == "canonical":
            transcript, rule = select_transcript(gene)
            if transcript is not None:
                exons = transcript.exons
                left, right = max(start, transcript.start), min(end, transcript.end)
        intervals = union_exons(exons)
        y = -lane
        visible_interval = [left, right] if left <= right else None
        if visible_interval is not None:
            ax.plot([left / 1e6, right / 1e6], [y, y], color=color, lw=0.65, zorder=2)
        visible_exons = []
        for a, b in intervals:
            a, b = max(start, a), min(end, b)
            if a <= b:
                visible_exons.append([a, b])
                ax.add_patch(
                    Rectangle(
                        (a / 1e6, y - 0.13),
                        (b - a + 1) / 1e6,
                        0.26,
                        facecolor=color,
                        edgecolor="none",
                        zorder=3,
                    )
                )
        if gene.strand in {"+", "-"} and right > left:
            locations = np.linspace(
                left, right, max(2, int((right - left) / (end - start + 1) * 18) + 1)
            )[1:-1]
            if not len(locations):
                locations = [(left + right) / 2]
            ax.scatter(
                np.asarray(locations) / 1e6,
                np.repeat(y, len(locations)),
                marker=">" if gene.strand == "+" else "<",
                s=9,
                color=color,
                linewidths=0,
                zorder=4,
            )
        label_center = (max(start, gene.start) + min(end, gene.end)) / 2e6
        label = ax.text(
            (left + right) / 2e6 if visible_interval is not None else label_center,
            y + 0.19,
            gene.name
            if visible_interval is not None
            else gene.name + " (selected transcript outside window)",
            fontsize=fontsize if target else fontsize - 0.5,
            weight="bold" if target else "normal",
            color=color,
            ha="center",
            va="bottom",
        )
        record = dict(
            gene_id=gene.gene_id,
            name=gene.name,
            target=target,
            lane=lane,
            lane_overlap_bp=scores[lane],
            color=color,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            visible_interval=visible_interval,
            union_exons=intervals,
            visible_exons=visible_exons,
            transcript=None if transcript is None else transcript.transcript_id,
            transcript_rule=rule,
            label_visible=True,
        )
        records.append(record)
        labels.append(label)
    used = max((r["lane"] for r in records), default=0) + 1
    ax.set_xlim(start / 1e6, end / 1e6)
    ax.set_ylim(-used + 0.45, 0.9)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelbottom=False, bottom=False)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    accepted = []
    for label, record in zip(labels, records):
        fits = _bound_text(label, ax, renderer)
        box = label.get_window_extent(renderer).expanded(1.05, 1.1)
        if not fits or any(box.overlaps(other) for other in accepted):
            if record["target"]:
                raise ValueError(
                    "Target gene label cannot fit; increase figure width or reduce fontsize"
                )
            label.set_visible(False)
            record["label_visible"] = False
        else:
            accepted.append(box)
    ax.gene_metadata = dict(
        genes=records,
        target_gene=target_gene,
        target_mode=target_mode,
        lanes=used,
        suppressed_labels=[r["gene_id"] for r in records if not r["label_visible"]],
        coordinate_convention="1-based inclusive",
        annotation=asdict(annotation),
    )
    return ax


def plot_ld_heatmap(
    ld,
    *,
    start,
    end,
    ax=None,
    max_snps=180,
    lead=None,
    colorbar=True,
    colorbar_ax=None,
    allowed_ids=None,
    cmap="white_red",
):
    """Draw a vector triangle with physical-coordinate cells and audited thinning.

    Only a real pairwise matrix can supply a heatmap. Repeated physical positions
    retain one deterministic SNP (lead preferred); rank thinning affects only
    this panel. Return the axis with ``ld_metadata`` and an optional colorbar.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    ld = as_ld(ld)
    if start < 1 or end <= start:
        raise ValueError("Heatmap requires a nonempty physical window")
    if ax is None:
        _, ax = plt.subplots(figsize=(7.09, 2.0))
    ax.set_xlim(start / 1e6, end / 1e6)
    ax.set_ylim((end - start) / 2e6, 0)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel("Position (Mb)")
    if not isinstance(ld, LDMatrix):
        reason = (
            "Pairwise LD unavailable"
            if ld is None
            else "Pairwise LD unavailable (lead vector only)"
        )
        ax.text(0.5, 0.5, reason, ha="center", va="center", transform=ax.transAxes, fontsize=7)
        ax.ld_metadata = dict(status="unavailable", reason=reason, selected_ids=[])
        return ax
    indices, metadata = heatmap_selection(ld, start, end, max_snps, lead, allowed_ids)
    if len(indices) < 2:
        ax.text(
            0.5,
            0.5,
            "Pairwise LD unavailable (<2 unique positions)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=7,
        )
        ax.ld_metadata = dict(metadata, status="unavailable")
        return ax
    pos = ld.positions[indices]
    edges = np.r_[
        max(start, pos[0] - (pos[1] - pos[0]) / 2),
        (pos[:-1] + pos[1:]) / 2,
        min(end, pos[-1] + (pos[-1] - pos[-2]) / 2),
    ]
    polygons, values = [], []
    for i, original_i in enumerate(indices):
        for j in range(i, len(indices)):
            a, b, c, d = edges[i], edges[i + 1], edges[j], edges[j + 1]
            vertices = [(a, c), (b, c), (b, d), (a, d)] if i != j else [(a, a), (b, b), (a, b)]
            polygons.append([((u + v) / 2e6, (v - u) / 2e6) for u, v in vertices])
            values.append(ld.r2[original_i, indices[j]])
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    if cmap == "white_red":
        cmap = LinearSegmentedColormap.from_list(
            "white_red", ["#ffffff", "#fee0d2", "#fc9272", "#de2d26", "#a50f15"]
        )
    norm = Normalize(0, 1)
    collection = PolyCollection(
        polygons,
        array=np.asarray(values),
        cmap=cmap,
        norm=norm,
        edgecolors="none",
        linewidths=0,
        rasterized=False,
        antialiased=False,
    )
    ax.add_collection(collection)
    if colorbar:
        bar = ax.figure.colorbar(
            collection,
            ax=ax if colorbar_ax is None else None,
            cax=colorbar_ax,
            ticks=LD_BINS,
            fraction=0.035,
            pad=0.02,
        )
        bar.set_label("r²", fontsize=7)
        if bar.solids is not None:
            bar.solids.set_rasterized(False)
        bar.ax.tick_params(labelsize=6.5, length=2)
        bar.outline.set_linewidth(0.4)
    ax.ld_metadata = dict(
        metadata,
        status="available",
        cell_edges_bp=edges.tolist(),
        selected_r2=ld.r2[np.ix_(indices, indices)].tolist(),
        vector=True,
        colormap=collection.cmap.name,
    )
    return ax


def _select_row(result, pair):
    if pair is not None:
        if not {"idx1", "idx2"} <= set(result.summary):
            raise ValueError("signal_pair only applies to a retained BF/SuSiE comparison")
        selected = result.summary[
            (result.summary.idx1 == pair[0]) & (result.summary.idx2 == pair[1])
        ]
        if len(selected) != 1:
            raise ValueError("signal_pair does not identify exactly one original component pair")
        return selected.iloc[0]
    if len(result.summary) > 1:
        raise ValueError("Multiple pairs require explicit original component IDs in signal_pair")
    return None if result.summary.empty else result.summary.iloc[0]


def _tracks(result, datasets, region, chromosome_policy):
    chromosome, start, end, build = region
    validate_region(chromosome, start, end, build)
    tracks = [_track(d) for d in datasets]
    coords = []
    for track in tracks:
        ids, pos, chrom, _, _ = track
        if normalize_chromosome(chrom, chromosome_policy) != normalize_chromosome(
            chromosome, chromosome_policy
        ):
            raise ValueError("Association chromosome disagrees with explicit region")
        if np.any(pos != np.floor(pos)) or np.any((pos < start) | (pos > end)):
            raise ValueError(
                "Every association coordinate must be an integer inside the explicit region"
            )
        if not set(result.results.get("snp", [])) <= set(ids):
            raise ValueError("Association dataset is missing analysed SNP IDs")
        coords.append(dict(zip(ids, pos)))
    for snp in set(coords[0]) & set(coords[1]):
        if coords[0][snp] != coords[1][snp]:
            raise ValueError("Shared SNP positions disagree")
    if "position" in result.results:
        for snp, pos in zip(result.results.snp, result.results.position):
            if coords[0][snp] != pos:
                raise ValueError("Association coordinates disagree with analysed result")
    return tracks


def _font():
    from matplotlib import font_manager

    names = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"):
        if name in names:
            return name
    return "DejaVu Sans"


def plot_coloc_locus(
    result,
    gwas,
    eqtl,
    *,
    region,
    annotation=None,
    target_gene=None,
    signal_pair=None,
    gwas_ld=None,
    eqtl_ld=None,
    heatmap_ld=None,
    leads=None,
    lead_rule="min_p",
    trait1="GWAS",
    trait2="eQTL",
    chromosome_policy="chr",
    output_prefix=None,
    formats=("svg", "pdf", "png"),
    figsize=(180 / 25.4, 6.6),
    fontsize=7.5,
    point_size=7,
    max_labels=2,
    max_heatmap_snps=180,
    target_mode="union",
    max_gene_lanes=6,
    heatmap_cmap="white_red",
    dpi=180,
):
    """Create a summary, two association tracks, genes and physical LD triangle.

    ``region=(chromosome, start, end, 'GRCh38')`` is explicit, 1-based inclusive.
    Association mappings/DataFrames require snp, chromosome, position and original
    pvalues (preferred), logp or z/beta/varbeta. Every supplied point is retained;
    all analysed IDs must be supplied. Counts distinguish displayed and analysed.
    LDMatrix/LeadLD or their field mappings are accepted separately for each
    association; heatmap_ld must be explicitly chosen, never borrowed silently.
    Leads default to minimum P (maximum plotted −log10P, ties position then ID);
    lead_rule='component' requires a selected row's hit1/hit2. Explicit two-item
    leads override this. Missing reference leads yield gray points, not zero LD.
    max_labels is 0–2 across both tracks; heatmap thinning never changes tracks.
    Multiple comparisons require original zero-based signal_pair IDs.

    Returns editable ``(figure, axes)`` with axes ordered summary/GWAS/eQTL/genes/LD.
    ``figure.coloc_metadata`` and the JSON sidecar contain full provenance and
    numeric selections. Output always includes a PNG preview; SVG text remains
    editable and PDF/heatmap stay vector. Only a local rc_context is changed.
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.text import Text

    if lead_rule not in {"min_p", "component"} or max_labels not in (0, 1, 2):
        raise ValueError("lead_rule is min_p/component; max_labels must be 0, 1 or 2")
    if any(f not in {"svg", "pdf", "png"} for f in formats):
        raise ValueError("Supported export formats are svg, pdf and png")
    row = _select_row(result, signal_pair)
    tracks = _tracks(result, (gwas, eqtl), region, chromosome_policy)
    chromosome, start, end, build = region
    chromosome = normalize_chromosome(chromosome, chromosome_policy)
    if annotation is not None:
        if (
            annotation.genome_build != build
            or normalize_chromosome(annotation.chromosome, chromosome_policy) != chromosome
        ):
            raise ValueError("Annotation build/chromosome mismatch")
        if (annotation.start, annotation.end) != (start, end):
            raise ValueError(
                "Prepare annotation for this exact window to retain complete overlap context"
            )
    elif target_gene is not None:
        raise ValueError("A target gene requires a prepared annotation")
    lds = [as_ld(gwas_ld), as_ld(eqtl_ld)]
    heatmap_ld = as_ld(heatmap_ld)
    for ld in (*lds, heatmap_ld):
        if ld is None:
            continue
        if "genome_build" in ld.provenance and ld.provenance["genome_build"] != build:
            raise ValueError("LD provenance genome build disagrees with region")
        if (
            "chromosome" in ld.provenance
            and normalize_chromosome(ld.provenance["chromosome"], chromosome_policy) != chromosome
        ):
            raise ValueError("LD provenance chromosome disagrees with region")
    if isinstance(heatmap_ld, LDMatrix):
        heat_coords = dict(zip(heatmap_ld.snp_ids, heatmap_ld.positions))
        for track in tracks:
            for snp, position in zip(track[0], track[1]):
                if snp in heat_coords and heat_coords[snp] != position:
                    raise ValueError("Heatmap and association coordinates disagree for " + snp)
    if leads is not None and len(leads) != 2:
        raise ValueError("leads must contain one SNP ID per association track")
    chosen = []
    for k, track in enumerate(tracks):
        ids, pos, _, y, _ = track
        if leads is not None:
            lead = str(leads[k])
        elif lead_rule == "component":
            lead = None if row is None else row.get(f"hit{k + 1}")
        else:
            lead = str(ids[np.lexsort((ids, pos, -y))[0]])
        if lead not in set(ids):
            raise ValueError(
                "Selected lead is absent from association data; choose an explicit lead or min_p rule"
            )
        chosen.append(lead)
    font = _font()
    rc = {
        "font.family": font,
        "font.size": fontsize,
        "axes.labelsize": fontsize,
        "xtick.labelsize": fontsize - 1,
        "ytick.labelsize": fontsize - 1,
        "axes.linewidth": 0.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "text.usetex": False,
    }
    with mpl.rc_context(rc):
        fig = plt.figure(figsize=figsize)
        grid = fig.add_gridspec(
            5,
            1,
            height_ratios=[0.86, 1.15, 1.15, 0.95, 1.55],
            left=0.10,
            right=0.91,
            top=0.975,
            bottom=0.075,
            hspace=0.15,
        )
        axes = [fig.add_subplot(grid[0])]
        axes.extend(
            fig.add_subplot(grid[j], sharex=None if j == 1 else axes[1]) for j in range(1, 5)
        )
        strip, ax1, ax2, genes_ax, ld_ax = axes
        strip.set_axis_off()
        posterior = None if row is None else {h: float(row[h]) for h in HYPOTHESES}
        eligible = posterior is not None and all(np.isfinite(v) for v in posterior.values())
        pair = None if row is None or "idx1" not in row else [int(row.idx1), int(row.idx2)]
        h4_status = None if row is None else row.get("h4_status")
        status = result.status if row is None else row.get("status", result.status)
        shared = len(set(tracks[0][0]) & set(tracks[1][0]))
        analysed = None if row is None else int(row.nsnps)
        weighted = None if row is None else any(w is not None for w in result.weights.values())
        target_name = (
            next((g.name for g in annotation.genes if matches_gene(g, target_gene)), target_gene)
            if annotation is not None
            else target_gene
        )
        title = f"{trait1} × {trait2}" + (f" | {target_name}" if target_name else "")
        strip.text(
            0, 1, title, va="top", fontsize=fontsize + 1, weight="bold", transform=strip.transAxes
        )
        location = f"{chromosome}:{start:,}–{end:,} {build} | {result.method}"
        location += (
            f" | n={analysed}"
            if analysed is not None
            else f" | coloc n unavailable; shared association n={shared}"
        )
        strip.text(0, 0.80, location, va="top", fontsize=fontsize - 0.5, transform=strip.transAxes)
        probabilities = (
            "   ".join(f"PP.H{j}={posterior[h]:.3g}" for j, h in enumerate(HYPOTHESES))
            if eligible
            else f"H0–H4 unavailable ({status})"
        )
        if h4_status == "posterior_underflow":
            probabilities += " (H4 underflow; conditional defined)"
        elif h4_status == "zero_shared_evidence":
            probabilities += " (zero shared evidence)"
        strip.text(0, 0.60, probabilities, va="top", fontsize=fontsize, transform=strip.transAxes)
        prior_text = "  ".join(
            f"{p}={result.priors[p]:g}" for p in ("p1", "p2", "p12") if p in result.priors
        )
        prior_text += f" | weighted priors: {'unavailable' if weighted is None else 'yes' if weighted else 'no'}"
        if pair is not None:
            prior_text += (
                f" | GWAS CS{pair[0] + 1} vs eQTL CS{pair[1] + 1}"
                if result.method == "coloc.susie"
                else f" | BF components {pair[0]}/{pair[1]}"
            )
        strip.text(0, 0.40, prior_text, va="top", fontsize=fontsize - 1, transform=strip.transAxes)
        handles = [
            Line2D([], [], marker="o", ls="", color=c, markersize=3.5, label=f"{a:g}–{b:g}")
            for c, a, b in zip(LD_COLORS, LD_BINS[:-1], LD_BINS[1:])
        ]
        handles += [
            Line2D(
                [], [], marker="o", ls="", color=MISSING_COLOR, markersize=3.5, label="LD missing"
            ),
            Line2D([], [], marker="D", ls="", color=LEAD_COLOR, markersize=3.5, label="Lead"),
        ]
        strip.legend(
            handles=handles,
            title="Lead r²",
            title_fontsize=fontsize - 1,
            fontsize=fontsize - 1,
            loc="lower left",
            bbox_to_anchor=(0, 0),
            frameon=False,
            ncol=7,
            borderaxespad=0,
            handletextpad=0.25,
            columnspacing=0.7,
            labelspacing=0.1,
        )
        track_meta = []
        cmap, norm = _scale()
        for k, (ax, track, ld, lead, name) in enumerate(
            zip((ax1, ax2), tracks, lds, chosen, (trait1, trait2))
        ):
            ids, pos, _, y, source = track
            r2, reason = lead_colors(ld, ids, pos, lead)
            colors = np.tile(np.asarray(mpl.colors.to_rgba(MISSING_COLOR)), (len(ids), 1))
            known = np.isfinite(r2)
            colors[known] = cmap(norm(r2[known]))
            ax.scatter(pos / 1e6, y, s=point_size, c=colors, linewidths=0, zorder=2)
            lead_i = int(np.flatnonzero(ids == lead)[0])
            ax.scatter(
                [pos[lead_i] / 1e6],
                [y[lead_i]],
                marker="D",
                s=point_size * 3.3,
                color=LEAD_COLOR,
                edgecolors="white",
                linewidths=0.4,
                zorder=4,
            )
            ax.set_ylim(0, max(1.0, float(np.max(y))) * 1.32)
            ax.set_yticks([tick for tick in ax.get_yticks() if 0 <= tick <= ax.get_ylim()[1]])
            ax.set_xlim(start / 1e6, end / 1e6)
            ax.set_ylabel("−log₁₀(P)")
            ax.tick_params(axis="x", labelbottom=False, bottom=False)
            ax.spines[["top", "right"]].set_visible(False)
            panel = "LD unavailable" if ld is None else ld.panel
            panel_title = ax.text(
                0.005,
                1.02,
                f"{name} · {panel} · n={len(ids)}",
                transform=ax.transAxes,
                va="bottom",
                fontsize=fontsize - 0.5,
            )
            if reason == "Selected lead absent from LD panel":
                ax.text(
                    1,
                    1.02,
                    "Lead absent from LD panel",
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=fontsize - 1,
                )
            label = None
            if k < max_labels:
                label = ax.annotate(
                    lead,
                    (pos[lead_i] / 1e6, y[lead_i]),
                    xytext=(pos[lead_i] / 1e6, y[lead_i] + ax.get_ylim()[1] * 0.10),
                    textcoords="data",
                    fontsize=fontsize - 1,
                    ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#555555"),
                    zorder=5,
                )
            record = dict(
                name=name,
                count=len(ids),
                snp_ids=ids.tolist(),
                positions=pos.astype(int).tolist(),
                minus_log10_p=y.tolist(),
                p_source=source,
                original_pvalues=None
                if "pvalues" not in (gwas, eqtl)[k]
                else np.asarray((gwas, eqtl)[k]["pvalues"]).tolist(),
                lead=lead,
                lead_rule="explicit" if leads is not None else lead_rule,
                panel=None if ld is None else ld.panel,
                ld_provenance=None if ld is None else ld.provenance,
                ld_kind=None if ld is None else ld.kind,
                r2=r2.tolist(),
                ld_status=reason,
                missing_ld_ids=ids[~known].tolist(),
                label_visible=label is not None,
            )
            track_meta.append(record)
            if label is not None:
                fig.canvas.draw()
                if not _bound_text(label, ax, fig.canvas.get_renderer()):
                    raise ValueError(
                        "Full SNP label does not fit; enlarge figure width or reduce fontsize"
                    )
                renderer = fig.canvas.get_renderer()
                title_box = panel_title.get_window_extent(renderer)
                label_box = Text.get_window_extent(label, renderer)
                if label_box.overlaps(title_box.expanded(1.02, 1.1)):
                    dx = title_box.x1 + 5 - label_box.x0
                    dy = 0
                    if label_box.x1 + dx > ax.bbox.x1 - 3:
                        dx = 0
                        dy = title_box.y0 - 4 - label_box.y1
                    old = ax.transData.transform(label.get_position())
                    label.set_position(ax.transData.inverted().transform(old + [dx, dy]))
                record["label_title_overlap"] = bool(
                    Text.get_window_extent(label, renderer).overlaps(title_box)
                )
        if annotation is None:
            genes_ax.set_axis_off()
            genes_ax.text(
                0.5,
                0.5,
                "Gene annotation unavailable (local annotation not supplied)",
                transform=genes_ax.transAxes,
                ha="center",
                fontsize=fontsize - 1,
            )
            genes_ax.gene_metadata = dict(status="unavailable")
        else:
            plot_gene_track(
                annotation,
                target_gene=target_gene,
                ax=genes_ax,
                target_mode=target_mode,
                max_lanes=max_gene_lanes,
                fontsize=fontsize - 0.5,
            )
        heatmap_lead = chosen[1]
        allowed_ids = set(result.results.get("snp", [])) or (set(tracks[0][0]) & set(tracks[1][0]))
        cax = (
            fig.add_axes(
                [0.929, ld_ax.get_position().y0 + 0.04, 0.012, ld_ax.get_position().height * 0.66]
            )
            if isinstance(heatmap_ld, LDMatrix)
            else None
        )
        plot_ld_heatmap(
            heatmap_ld,
            start=start,
            end=end,
            ax=ld_ax,
            max_snps=max_heatmap_snps,
            lead=heatmap_lead,
            colorbar_ax=cax,
            allowed_ids=allowed_ids,
            cmap=heatmap_cmap,
        )
        if cax is not None and ld_ax.ld_metadata["status"] != "available":
            cax.remove()
        ld_ax.set_xlabel(f"{chromosome} position (Mb) · {build}")
        from matplotlib.ticker import MaxNLocator

        ticks = MaxNLocator(nbins=6).tick_values(start / 1e6, end / 1e6)
        ld_ax.set_xticks(ticks[(ticks >= start / 1e6) & (ticks <= end / 1e6)])
        if heatmap_ld is not None:
            ld_ax.text(
                0.005,
                0.05,
                f"Pairwise LD · {heatmap_ld.panel}",
                transform=ld_ax.transAxes,
                fontsize=fontsize - 1,
                va="bottom",
                bbox=dict(facecolor="white", edgecolor="none", pad=1),
            )
        conditional = None if row is None else row.get("snp_pp_column", "SNP.PP.H4")
        candidate = None
        if conditional in result.results and np.isfinite(result.results[conditional]).any():
            best = result.results.sort_values([conditional, "snp"], ascending=[False, True]).iloc[0]
            candidate = dict(
                snp=str(best.snp), conditional_posterior=float(best[conditional]), field=conditional
            )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        extents = []
        for text in fig.findobj(Text):
            if not text.get_visible() or not text.get_text():
                continue
            box = text.get_window_extent(renderer)
            inside = (
                box.x0 >= -0.5
                and box.y0 >= -0.5
                and box.x1 <= fig.bbox.width + 0.5
                and box.y1 <= fig.bbox.height + 0.5
            )
            extents.append(
                dict(text=text.get_text(), bounds=list(box.bounds), inside_canvas=inside)
            )
        metadata = dict(
            schema="pycoloc-locus-1",
            method=result.method,
            status=status,
            h4_status=h4_status,
            signal_pair=pair,
            component_indices="original zero-based; CS label = original index + 1",
            posteriors=posterior,
            priors=result.priors,
            weighted_priors=weighted,
            analysed_snps=analysed,
            shared_association_snps=shared,
            region=dict(chromosome=chromosome, start=start, end=end, genome_build=build),
            tracks=track_meta,
            genes=genes_ax.gene_metadata,
            heatmap=ld_ax.ld_metadata,
            shared_candidate=candidate,
            font=font,
            font_size=fontsize,
            figsize=list(figsize),
            rasterized_heatmap=False,
            renderer_extents=extents,
            extent_pass=all(t["inside_canvas"] for t in extents),
        )
        fig.coloc_metadata = _clean(metadata)
        if output_prefix is not None:
            if not fig.coloc_metadata["extent_pass"]:
                outside = [t["text"] for t in extents if not t["inside_canvas"]]
                plt.close(fig)
                raise ValueError(
                    "Layout extends outside figure; enlarge figsize or shorten labels: "
                    + repr(outside)
                )
            prefix = Path(output_prefix)
            prefix.parent.mkdir(parents=True, exist_ok=True)
            outputs = []
            for extension in dict.fromkeys((*formats, "png")):
                path = Path(str(prefix) + "." + extension)
                fig.savefig(path, dpi=dpi, facecolor="white")
                outputs.append(str(path.resolve()))
            fig.coloc_metadata["outputs"] = outputs
            Path(str(prefix) + ".json").write_text(
                json.dumps(fig.coloc_metadata, indent=2, allow_nan=False) + "\n"
            )
        return fig, axes
