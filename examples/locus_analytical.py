"""Small analytical plot illustration; not a real locus or simulated study."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np

from pycoloc import LDMatrix, coloc_abf, plot_coloc_locus
from pycoloc.annotation import ExonModel, GeneModel, RegionAnnotation

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-prefix", type=Path, required=True)
args = parser.parse_args()
ids = ["a", "b", "c"]
positions = [100, 200, 300]
trait1 = dict(
    snp=ids,
    position=positions,
    chromosome="22",
    beta=[0.1, 0.5, -0.2],
    varbeta=[0.01, 0.02, 0.03],
    type="cc",
)
trait2 = dict(trait1, beta=[0.0, 0.4, 0.1], type="quant", sdY=1.2)
result = coloc_abf(trait1, trait2)
annotation = RegionAnnotation(
    "chr22",
    50,
    350,
    "GRCh38",
    (
        GeneModel(
            "illustration",
            "Illustrative gene",
            "chr22",
            80,
            320,
            "+",
            exons=(
                ExonModel("exon1", "chr22", 90, 140, "+"),
                ExonModel("exon2", "chr22", 240, 300, "+"),
            ),
        ),
    ),
    {"source": "hand-written analytical example; not biological annotation"},
)
ld = LDMatrix(
    np.array([[1.0, 0.2, 0.0], [0.2, 1.0, 0.1], [0.0, 0.1, 1.0]]),
    ids,
    positions,
    kind="r",
    panel="Analytical matrix",
    provenance={"source": "hand-written"},
)
fig, axes = plot_coloc_locus(
    result,
    trait1,
    trait2,
    region=("22", 50, 350, "GRCh38"),
    annotation=annotation,
    target_gene="illustration",
    gwas_ld=ld,
    eqtl_ld=ld,
    heatmap_ld=ld,
    trait1="Analytical trait 1",
    trait2="Analytical trait 2",
    output_prefix=args.output_prefix,
    max_labels=0,
)
assert len(axes) == 5
for suffix in [".svg", ".pdf", ".png", ".json"]:
    assert args.output_prefix.with_suffix(suffix).is_file()
print("Saved analytical illustration:", args.output_prefix)
