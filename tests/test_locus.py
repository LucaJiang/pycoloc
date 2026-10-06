"""Analytic LD, geometry, identity and export tests for the public feature."""

import json
from dataclasses import replace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pycoloc import (
    ColocResult,
    LDMatrix,
    LeadLD,
    coloc_abf,
    coloc_susie,
    plot_coloc_locus,
    plot_gene_track,
    plot_ld_heatmap,
)
from pycoloc.annotation import ExonModel, GeneModel, RegionAnnotation, TranscriptModel
from pycoloc.annotation.models import select_transcript


@pytest.fixture
def annotation():
    e1, e2 = ExonModel("e1", "chr22", 50, 150, "+"), ExonModel("e2", "chr22", 220, 260, "+")
    e3 = ExonModel("e3", "chr22", 145, 220, "+")
    t1 = TranscriptModel("long", "chr22", 50, 500, "+", (e1, e2))
    t2 = TranscriptModel("mane", "chr22", 145, 400, "+", (e3,), ("MANE_Select",))
    g = GeneModel("g.1", "Target", "chr22", 50, 600, "+", (t1, t2))
    context = GeneModel(
        "g2", "Context", "chr22", 350, 500, "-", (), (ExonModel("e4", "chr22", 380, 400, "-"),)
    )
    return RegionAnnotation("chr22", 90, 450, "GRCh38", (g, context), {"source": "test"})


def test_gene_union_canonical_highlight(annotation):
    ax = plot_gene_track(annotation, target_gene="g")
    target = ax.gene_metadata["genes"][0]
    assert target["target"] and target["name"] == "Target" and target["color"] == "#222222"
    assert target["union_exons"] == ((50, 260),) and target["visible_exons"] == [[90, 260]]
    assert target["start"] == 50 and target["end"] == 600
    assert len(ax.gene_metadata["genes"]) == 2
    plt.close(ax.figure)
    ax = plot_gene_track(annotation, target_gene="g", target_mode="canonical")
    assert ax.gene_metadata["genes"][0]["transcript"] == "mane"
    assert ax.gene_metadata["genes"][0]["transcript_rule"] == "annotated MANE_Select"
    plt.close(ax.figure)
    g = annotation.genes[0]
    representative, rule = select_transcript(
        replace(g, transcripts=tuple(replace(t, tags=()) for t in g.transcripts))
    )
    assert representative.transcript_id == "long" and rule.startswith("representative fallback")


def test_signed_r_geometry_and_permutation():
    r = np.array([[1, -0.5, 0], [-0.5, 1, 0.8], [0, 0.8, 1.0]])
    ids = np.array(["a", "b", "c"])
    pos = np.array([100, 300, 900])
    ld = LDMatrix(r, ids, pos, "r", "panel")
    ax = plot_ld_heatmap(ld, start=50, end=1000)
    np.testing.assert_allclose(ax.collections[0].get_array(), [1, 0.25, 0, 1, 0.64, 1])
    np.testing.assert_allclose(
        ax.collections[0].get_paths()[1].vertices[:4] * 1e6,
        [[125, 75], [200, 0], [400, 200], [325, 275]],
    )
    assert ax.ld_metadata["cell_edges_bp"] == [50, 200, 600, 1000]
    assert len(ax.figure.axes) == 2 and ax.figure.axes[1].get_ylabel() == "r²"
    assert not ax.collections[0].get_rasterized()
    plt.close(ax.figure)
    permutation = [2, 0, 1]
    shuffled = LDMatrix(
        (r * r)[np.ix_(permutation, permutation)], ids[permutation], pos[permutation], "r2", "panel"
    )
    ax = plot_ld_heatmap(shuffled, start=50, end=1000)
    assert ax.ld_metadata["selected_ids"] == ["a", "b", "c"]
    assert ax.ld_metadata["original_indices"] == [1, 2, 0]
    np.testing.assert_allclose(ax.collections[0].get_array(), [1, 0.25, 0, 1, 0.64, 1])
    plt.close(ax.figure)


def test_thinning_lead_and_repeated_positions():
    ids = np.array(list("abcdefghij"))
    positions = np.array([100, 200, 200, 300, 400, 500, 600, 700, 800, 900])
    ld = LDMatrix(np.eye(10), ids, positions, "r2", "test")
    ax = plot_ld_heatmap(ld, start=50, end=950, max_snps=3, lead="c")
    assert ax.ld_metadata["selected_ids"] == ["c", "f", "j"]
    assert ax.ld_metadata["repeated_position_omissions"] == ["b"]
    assert ax.ld_metadata["positions"] == [200, 500, 900] and ax.ld_metadata["lead_retained"]
    plt.close(ax.figure)


@pytest.mark.parametrize(
    "change,match",
    [
        ({"kind": "unspecified"}, "kind"),
        ({"values": [[1, 0.1], [0.2, 1]]}, "symmetric"),
        ({"values": [[0.9, 0], [0, 1]]}, "diagonal"),
        ({"values": [[1, 1.2], [1.2, 1]]}, "range"),
        ({"values": [[1, np.nan], [np.nan, 1]]}, "finite"),
        ({"snp_ids": ["a", "a"]}, "unique"),
        ({"positions": [1, 1.2]}, "integer"),
        ({"values": [[1]]}, "shape"),
    ],
)
def test_invalid_ld(change, match):
    kw = dict(values=np.eye(2), snp_ids=["a", "b"], positions=[100, 200], kind="r", panel="test")
    with pytest.raises(ValueError, match=match):
        LDMatrix(**{**kw, **change})


def test_full_tracks_missing_per_panel_and_exports(datasets, annotation, tmp_path):
    d1, d2 = datasets
    result = coloc_abf(d1, d2, prior_weights1=[1, 9, 2, 1])
    ld1 = LDMatrix([[1, -0.5], [-0.5, 1]], ["b", "c"], [200, 300], "r", "FIN")
    ld2 = LDMatrix([[1, 0.8], [0.8, 1]], ["b", "c"], [200, 300], "r2", "EUR")
    prefix = tmp_path / "locus"
    before = matplotlib.rcParams.copy()
    fig, axes = plot_coloc_locus(
        result,
        pd.DataFrame({k: v for k, v in d1.items() if k not in ["type"]}),
        d2,
        region=("22", 90, 450, "GRCh38"),
        annotation=annotation,
        target_gene="g",
        gwas_ld=ld1,
        eqtl_ld=ld2,
        heatmap_ld=ld2,
        output_prefix=prefix,
    )
    assert matplotlib.rcParams == before
    assert len(axes) == 5
    meta = fig.coloc_metadata
    assert meta["tracks"][0]["r2"] == [None, 1, 0.25, None]
    assert meta["tracks"][1]["r2"] == [None, 1, 0.8, None]
    assert meta["tracks"][0]["missing_ld_ids"] == ["a", "d"]
    assert (
        meta["shared_candidate"]["snp"]
        == result.results.loc[result.results["SNP.PP.H4"].idxmax(), "snp"]
    )
    for ax, d in zip(axes[1:3], datasets):
        np.testing.assert_array_equal(
            ax.collections[0].get_offsets()[:, 0], np.array(d["position"]) / 1e6
        )
        np.testing.assert_allclose(
            ax.collections[0].get_offsets()[:, 1], -np.log10(d["pvalues"]), atol=0, rtol=0
        )
    assert meta["extent_pass"]
    assert all(
        prefix.with_suffix("." + ext).stat().st_size > 100 for ext in ["svg", "pdf", "png", "json"]
    )
    svg = prefix.with_suffix(".svg").read_text()
    assert "<text" in svg and "<image" not in svg
    assert json.loads(prefix.with_suffix(".json").read_text())["genes"]["genes"][0]["target"]
    plt.close(fig)


def test_unknown_lead_missing_vector_and_no_cs(datasets):
    d1, d2 = datasets
    result = ColocResult(
        pd.DataFrame(),
        pd.DataFrame({"snp": []}),
        {"p1": 1e-4, "p2": 1e-4, "p12": 5e-6},
        "coloc.susie",
        "no_cs",
    )
    matrix = LDMatrix(np.eye(2), ["c", "d"], [300, 400], "r", "FIN")
    vector = LeadLD([1, 0.3], ["b", "c"], "b", "r2", "EUR")
    fig, axes = plot_coloc_locus(
        result,
        d1,
        d2,
        region=("22", 90, 450, "GRCh38"),
        gwas_ld=matrix,
        eqtl_ld=vector,
        heatmap_ld=vector,
    )
    meta = fig.coloc_metadata
    assert meta["analysed_snps"] is None and meta["shared_association_snps"] == 4
    assert meta["posteriors"] is None and meta["signal_pair"] is None
    assert meta["tracks"][0]["r2"] == [None] * 4 and "absent" in meta["tracks"][0]["ld_status"]
    assert meta["heatmap"]["status"] == "unavailable"
    assert any("unavailable (no_cs)" in t.get_text() for t in axes[0].texts)
    plt.close(fig)
    with pytest.raises(ValueError, match="reference differs"):
        plot_coloc_locus(
            result, d1, d2, region=("22", 90, 450, "GRCh38"), eqtl_ld=vector, leads=("b", "c")
        )


def test_original_components_and_required_pair(datasets, fits):
    result = coloc_susie(*fits)
    with pytest.raises(ValueError, match="Multiple pairs"):
        plot_coloc_locus(result, *datasets, region=("22", 90, 450, "GRCh38"))
    fig, axes = plot_coloc_locus(
        result,
        *datasets,
        region=("22", 90, 450, "GRCh38"),
        signal_pair=(2, 1),
        lead_rule="component",
    )
    assert fig.coloc_metadata["signal_pair"] == [2, 1]
    assert any("GWAS CS3 vs eQTL CS2" in t.get_text() for t in axes[0].texts)
    assert fig.coloc_metadata["tracks"][0]["lead"] == "c"
    plt.close(fig)


def test_coordinate_conflicts_rejected(datasets):
    result = coloc_abf(*datasets)
    with pytest.raises(ValueError, match="Shared SNP positions"):
        plot_coloc_locus(
            result,
            datasets[0],
            {**datasets[1], "position": [101, 200, 300, 400]},
            region=("22", 90, 450, "GRCh38"),
        )
    ld = LDMatrix(np.eye(2), ["b", "c"], [201, 300], "r", "FIN")
    with pytest.raises(ValueError, match="coordinates disagree"):
        plot_coloc_locus(result, *datasets, region=("22", 90, 450, "GRCh38"), gwas_ld=ld)


@pytest.mark.parametrize(
    "start,end,strand",
    [
        (10.5, 20, "+"),
        (float("nan"), 20, "+"),
        (10, float("inf"), "+"),
        (20, 10, "+"),
        (0, 10, "+"),
        (True, 20, "+"),
        (10, 20, "invalid"),
    ],
)
def test_invalid_annotation_models(start, end, strand):
    with pytest.raises(ValueError):
        ExonModel("e", "chr22", start, end, strand)


def test_annotation_parent_chromosome_and_coordinate_boundaries():
    e = ExonModel("e", "chr21", 100, 150, "+")
    with pytest.raises(ValueError, match="chromosome"):
        TranscriptModel("t", "chr22", 50, 200, "+", (e,))
    e = ExonModel("e", "chr22", 100, 250, "+")
    with pytest.raises(ValueError, match="outside parent"):
        TranscriptModel("t", "chr22", 50, 200, "+", (e,))


def test_canonical_transcript_outside_window():
    inside = TranscriptModel("in", "chr22", 50, 200, "+", (ExonModel("a", "chr22", 60, 100, "+"),))
    outside = TranscriptModel(
        "out", "chr22", 500, 600, "+", (ExonModel("b", "chr22", 500, 550, "+"),), ("MANE_Select",)
    )
    gene = GeneModel("g", "Target", "chr22", 50, 600, "+", (inside, outside))
    a = RegionAnnotation("chr22", 90, 450, "GRCh38", (gene,))
    ax = plot_gene_track(a, target_gene="g", target_mode="canonical")
    assert len(ax.lines) == 0 and len(ax.patches) == 0
    assert ax.gene_metadata["genes"][0]["visible_interval"] is None
    assert ax.gene_metadata["genes"][0]["transcript"] == "out"
    assert "outside window" in ax.texts[0].get_text()
    plt.close(ax.figure)


def test_ld_id_inputs_are_copied():
    ids = np.array(["a", "b"])
    ld = LDMatrix(np.eye(2), ids, [100, 200], "r", "test")
    vector = LeadLD([1, 0.2], ids, "a", "r2", "test")
    ids[:] = ["b", "a"]
    assert list(ld.snp_ids) == ["a", "b"] and list(vector.snp_ids) == ["a", "b"]


def test_heatmap_shared_support_original_indices_and_position_validation(datasets):
    result = coloc_abf(*datasets)
    # Extra reference SNP is in the region but absent from both analysed tracks.
    ld = LDMatrix(np.eye(5), ["extra", "d", "b", "a", "c"], [250, 400, 200, 100, 300], "r", "EUR")
    fig, axes = plot_coloc_locus(
        result, *datasets, region=("22", 90, 450, "GRCh38"), heatmap_ld=ld, max_heatmap_snps=4
    )
    assert fig.coloc_metadata["heatmap"]["selected_ids"] == ["a", "b", "c", "d"]
    assert fig.coloc_metadata["heatmap"]["original_indices"] == [3, 2, 4, 1]
    assert fig.coloc_metadata["heatmap"]["omitted_not_shared"] == ["extra"]
    assert len(axes[1].collections[0].get_offsets()) == 4
    plt.close(fig)
    bad = LDMatrix(np.eye(2), ["a", "b"], [101, 200], "r", "EUR")
    with pytest.raises(ValueError, match="Heatmap and association coordinates"):
        plot_coloc_locus(result, *datasets, region=("22", 90, 450, "GRCh38"), heatmap_ld=bad)


@pytest.mark.parametrize("provenance", [{"genome_build": "GRCh37"}, {"chromosome": "chr21"}])
def test_ld_provenance_mismatch(datasets, provenance):
    ld = LDMatrix(np.eye(2), ["a", "b"], [100, 200], "r", "test", provenance)
    with pytest.raises(ValueError, match="LD provenance"):
        plot_coloc_locus(
            coloc_abf(*datasets), *datasets, region=("22", 90, 450, "GRCh38"), heatmap_ld=ld
        )


def test_clipped_export_is_rejected(datasets, tmp_path):
    prefix = tmp_path / "must_not_publish"
    with pytest.raises(ValueError, match="Layout extends outside"):
        plot_coloc_locus(
            coloc_abf(*datasets),
            *datasets,
            region=("22", 90, 450, "GRCh38"),
            trait1="A deliberately unbounded title " * 80,
            output_prefix=prefix,
        )
    assert not list(tmp_path.glob("must_not_publish.*"))
