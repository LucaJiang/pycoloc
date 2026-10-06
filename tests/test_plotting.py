import copy

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycoloc import coloc_abf, coloc_susie, plot_coloc


def test_shared_candidate_uses_corrected_main_column(datasets):
    result = coloc_abf(*datasets, prior_weights1=[1e8, 1, 1, 1])
    # The BF-only lead is b, but the independently supplied prior selects a.
    assert result.results.loc[result.results["internal.sum.lABF"].idxmax(), "snp"] == "b"
    fig, axes = plot_coloc(result, *datasets)
    assert fig.coloc_metadata["shared_candidate"] == "a"
    assert fig.coloc_metadata["shared_posterior_column"] == "SNP.PP.H4"
    for ax in axes:
        shared = [c for c in ax.collections if c.get_label() == "Shared candidate | H4"]
        assert len(shared) == 1 and shared[0].get_offsets()[0, 0] == 0.0001
    plt.close(fig)


def test_title_probabilities_coordinates_and_original_p(datasets):
    result = coloc_abf(*datasets, prior_weights1=[1, 2, 3, 4])
    before = matplotlib.rcParams.copy()
    fig, axes = plot_coloc(result, *datasets, trait1="GWAS", trait2="Gene")
    for i in range(5):
        assert (
            f"PP.H{i} = {result.summary.iloc[0][f'PP.H{i}.abf']:.6g}"
            in fig.coloc_metadata["probabilities"]
        )
    for ax, d in zip(axes, datasets):
        assert_allclose(ax.collections[0].get_offsets()[:, 0], np.array(d["position"]) / 1e6)
        assert_allclose(ax.collections[0].get_offsets()[:, 1], -np.log10(d["pvalues"]))
    assert any("variant-specific priors: yes" in t.get_text() for t in fig.texts)
    assert fig.coloc_metadata["p_sources"] == ["original summary P"] * 2
    assert before == matplotlib.rcParams
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for t in fig.texts:
        box = t.get_window_extent(renderer)
        assert (
            box.x0 >= 0 and box.x1 <= fig.bbox.width and box.y0 >= 0 and box.y1 <= fig.bbox.height
        )
    plt.close(fig)


def test_pair_selection_and_cs(fits, datasets):
    result = coloc_susie(*fits)
    with pytest.raises(ValueError, match="multiple signal pairs"):
        plot_coloc(result, *datasets)
    fig, axes = plot_coloc(result, *datasets, signal_pair=(2, 1))
    assert fig.coloc_metadata["signal_pair"] == (2, 1)
    assert fig.coloc_metadata["leads"] == ["c", "c"]
    assert "components (2, 1)" in fig.texts[0].get_text()
    assert axes[0].collections[1].get_offsets()[0, 0] == 0.0003
    plt.close(fig)


def test_empty_and_nonconverged_are_explicit(fits, datasets):
    f1, f2 = copy.deepcopy(fits)
    f1["sets"] = dict(cs=[], cs_index=[])
    f1["converged"] = False
    with pytest.warns(UserWarning):
        result = coloc_susie(f1, f2)
    fig, _ = plot_coloc(result, *datasets)
    text = " ".join(t.get_text() for t in fig.texts)
    assert "no_cs" in text and "not converged" in text and "PP.H4 = 0" not in text
    plt.close(fig)


def test_missing_data_rejected(datasets):
    result = coloc_abf(*datasets)
    d = copy.deepcopy(datasets[0])
    del d["position"]
    with pytest.raises(ValueError, match="position"):
        plot_coloc(result, d, datasets[1])
    with pytest.raises(ValueError, match="same chromosome"):
        plot_coloc(result, dict(datasets[0], chromosome="21"), datasets[1])
    with pytest.raises(ValueError, match="conflicting plot positions"):
        plot_coloc(result, dict(datasets[0], position=[101, 200, 300, 400]), datasets[1])
    with pytest.raises(ValueError, match="analyzed positions"):
        plot_coloc(
            result,
            dict(datasets[0], position=[101, 200, 300, 400]),
            dict(datasets[1], position=[101, 200, 300, 400]),
        )


def test_zero_p_recovery_and_clipping(datasets):
    result = coloc_abf(*datasets)
    d = copy.deepcopy(datasets[0])
    d["pvalues"][0] = 0
    fig, ax = plot_coloc(result, d, datasets[1])
    assert "recovered" in fig.coloc_metadata["p_sources"][0]
    assert np.isfinite(ax[0].collections[0].get_offsets()).all()
    plt.close(fig)
    del d["beta"], d["varbeta"]
    fig, ax = plot_coloc(result, d, datasets[1])
    assert "clipped" in fig.coloc_metadata["p_sources"][0]
    assert ax[0].collections[0].get_offsets()[0, 1] > 300
    plt.close(fig)
