# Offline, publication-size locus plots (0.2.0)

`plot_coloc_locus` returns an editable matplotlib figure and five axes: compact
coloc summary, GWAS association, eQTL association, genes and physical-coordinate
pairwise LD. `plot_coloc` remains available unchanged. The new plotting API does
not alter `ColocResult` schema 2 or the corrected 0.1.1 statistical calculations.

```python
from pycoloc import LDMatrix, prepare_region_annotation, plot_coloc_locus

# All input objects and metadata are supplied by the caller.
annotation = prepare_region_annotation(
    local_gtf, chromosome=chromosome, start=start, end=end,
    genome_build="GRCh38", version=annotation_version,
    cache_dir="annotation-index", source_name="GENCODE",
)
ld1 = LDMatrix(r1, ids1, positions1, kind="r", panel="Trait 1 reference",
               provenance=provenance1)
ld2 = LDMatrix(r2, ids2, positions2, kind="r", panel="Trait 2 reference",
               provenance=provenance2)
fig, axes = plot_coloc_locus(
    result, dataset1, dataset2, region=(chromosome, start, end, "GRCh38"),
    annotation=annotation, target_gene=target_gene, signal_pair=signal_pair,
    gwas_ld=ld1, eqtl_ld=ld2, heatmap_ld=ld2, output_prefix="figures/locus",
)
```

Supply your existing result, association data and actual LD in this example;
pycoloc does not compute LD or invent allele coding. Mappings and DataFrames
require `snp`, `chromosome`, integer `position` and `pvalues` (preferred), `logp`,
`z` or `beta`/`varbeta`. Coordinates are positive **1-based inclusive**, including
GTF/GFF3 models. Every supplied association point must lie in the explicit window;
every analysed SNP must be present. Shared IDs must have matching positions in
both traits, the result, and supplied LD. The plot displays all supplied points,
which can exceed the common analysed set; both counts are recorded. It never
filters associations to LD support. Original P values are retained in metadata;
P=0 recovery from logP/Wald statistics or display clipping is explicitly identified.

The summary uses one exact comparison, with all five posteriors and actual priors.
Multiple pairs require explicit original zero-based component IDs. CS labels use
original component index + 1, not a re-numbered retained-set ordinal. No-CS/empty
results show unavailable posteriors and analysed count; shared association count
remains available. True zero shared evidence and floating H4 underflow are labeled
separately. The corrected main conditional SNP posterior determines the recorded
shared candidate; it is never confused with the joint H4/SNP probability.

Default leads are minimum original P (or maximum displayed −log10P when derived),
with position then SNP-ID tie breaks. `lead_rule="component"` uses selected
`hit1`/`hit2`; `leads=(gwas_id, eqtl_id)` overrides the rule explicitly. The figures
remain **marginal** association tracks, not residual/conditional SuSiE association
plots. A reference lead absent from its LD panel produces gray colors and an
explicit caption, without selecting a replacement silently.

## Annotation source, cache and fallback policy

Install the optional local indexer with `python -m pip install '.[annotation]'` from this repository.
Ordinary colocalisation and rendering a prepared `RegionAnnotation` require no
gffutils, R, prusie or network access. `GeneModel`, `TranscriptModel`, `ExonModel`
and `RegionAnnotation` are source/matplotlib-independent dataclasses in
`pycoloc.annotation`; `annotation_to_dict`/`annotation_from_dict` save/load full
prepared models offline.

Human GRCh38 local GENCODE comprehensive GTF is the default supported source.
Acquire a **specific release explicitly** from [GENCODE Human](https://www.gencodegenes.org/human/)
and retain its source URL, release and SHA256. Plotting never downloads data.
There is no remote Ensembl fallback: a missing/corrupt local file is an explicit
error or cache rebuild from the supplied source, never a different annotation.
Other complete local GTF/GFF3 files can be indexed with an explicit source/version;
explicit gene and transcript/mRNA records are required (inference is disabled).
The build must be declared GRCh38; a conflicting source header is rejected.

The SQLite backend uses [gffutils create_db](https://daler.github.io/gffutils/autodocs/gffutils.create.create_db.html)
with inference disabled, explicit duplicate errors and actual GTF/GFF3 parent
relationships. [Overlap queries](https://daler.github.io/gffutils/autodocs/gffutils.interface.FeatureDB.region.html)
include genes crossing the window boundaries and retain all complete transcript
and exon coordinates. Only rendering clips the visible segments. Default exons
are the union across a gene's transcripts, merging overlapping/adjacent inclusive
intervals. `target_mode="canonical"` prioritizes MANE_Select, then
Ensembl_canonical; otherwise it records a representative fallback chosen by
longest union-exon length and transcript ID. It never calls an untagged transcript
annotated canonical. Context genes remain collapsed unions, with directional marks
on both strands. Target identity can be a gene ID (with or without version) or
an unambiguous gene name.

Cache identity includes indexed-file SHA256/path, source/version/build, indexer
version, schema, chromosome policy, scope and upstream provenance. Valid cache
metadata, database SHA256 and SQLite integrity must agree before reuse. A changed
identity creates a separate index; corrupted cache files are rebuilt. POSIX file
locks serialize same-index builds and temporary files plus atomic replacements
prevent partial indices from being reused. The first index can take time and disk
space for a whole genome; no performance claim is made.

Chromosome spelling normalization is explicit: `chromosome_policy="chr"` adds a
leading `chr`, `"none"` removes it, and `"exact"` requires exact spelling. Ambiguous
normalizations fail. For a deliberate whole-chromosome subset, set
`index_scope=["chr22"]` and preserve the full-source hash and deterministic subset
rule in `source_provenance`. Requests outside that scope fail. Coordinate-window
subsets are not supported because boundary-spanning models would be incomplete.
A chromosome subset is not a genome-wide cache. Preserve the original full-source identity and the subset rule.

## LD, geometry and missing information

`LDMatrix(values, snp_ids, positions, kind, panel, provenance={})` validates unique
IDs, square dimensions, finite range, symmetry and unit diagonal at absolute
roundoff tolerance 1e-8 (configurable up to 1e-4). `kind` is explicitly `"r"` or
`"r2"`; signed r is squared exactly once. Metadata retains kind and provenance.
Use separate FIN and EUR objects for corresponding tracks. Select the single
heatmap panel explicitly with `heatmap_ld`; it is not silently pooled or borrowed.
Sample counts identify founders, not independently verified unrelated subjects.

`LeadLD(values, snp_ids, lead_id, kind, panel)` can color associations but never
creates a pairwise heatmap. Its reference lead must match the requested lead.
Missing LD is gray, separate from the blue 0–0.2 bin. Discrete boundaries are
0, .2, .4, .6, .8, 1; bins are left-inclusive and the final bin includes 1.
The highlighted purple diamond is independent of the LD color scale.

Only the heatmap can be thinned. SNPs are sorted by position then ID, repeated
positions keep a deterministic representative (requested lead preferred), and
at most `max_heatmap_snps` position-ranked evenly spaced SNPs remain, preserving
the heatmap lead (the eQTL lead in the composite plot). All selected IDs, exact
positions, original matrix indices, omissions and rules are recorded. Matrix
rows/columns are reindexed together; no association order is changed. Physical
cell edges are midpoint boundaries between selected positions, transformed by
x=(u+v)/2, y=(v−u)/2. This is a physical genomic triangle, not equal SNP-index
spacing; thinning necessarily coarsens displayed cells. Repeated-position
omissions are explicitly distinguishable from rank thinning and out-of-window LD.

## Styling, export and local example bundles

Default width is 180 mm, title 8.5 pt, axes 7.5 pt and secondary labels 6.5–7 pt.
The five-panel grid replaces the historical large multi-line header/two scatter
tracks with a compact summary/legend, identity-bound per-trait LD colors, complete
gene context and an audited vector triangle. Arial/Helvetica are preferred when
installed, then Liberation Sans/DejaVu Sans; metadata records the actual font.
Style changes use local `rc_context`, never global matplotlib settings.

Useful knobs are `figsize`, `fontsize`, `point_size`, `max_labels` (0–2 total),
`max_heatmap_snps`, `max_gene_lanes`, `target_mode` and explicit `leads`.
Context gene labels that collide are suppressed and recorded; gene models stay.
The target is prioritized. Lead labels are bounded and avoid their panel captions.
Especially dense loci may still need manual label placement in the editable SVG
or a wider figure. Metadata reports renderer bounds; inspect the PNG at intended
print size before publication. `plot_gene_track(annotation, ...)` and
`plot_ld_heatmap(ld, start=..., end=...)` also work independently, returning axes
with `gene_metadata` and `ld_metadata` respectively.

`output_prefix` saves editable SVG, vector PDF, PNG preview and a JSON sidecar;
PNG is always produced even if `formats` requests only vector formats. Heatmap
polygons are not rasterized. JSON includes full annotation provenance/models,
P values/coordinates, tracks, exact component pair, corrected candidate,
posteriors/priors, LD order/subset, labels and font/renderer bounds.

A user-supplied local bundle can reproduce figures with the installed package:

```sh
python examples/locus_region.py --config example_bundle.json \
  --output-prefix figures/locus --max-heatmap-snps 140
```

The bundle is JSON with `method` (`abf` or `susie`), `datasets` (two JSON mappings),
`region` (`chromosome`, `start`, `end`, `genome_build`), `target_gene`, optional
`coloc_options`, `signal_pair`, `leads`/`lead_rule`, `trait1`/`trait2` and
`plot_options`. SuSiE bundles include two `fits` entries (`json` metadata and
`npz` arrays); no refitting is done. `annotation_snapshot` points to a prepared
model, or `annotation_source` supplies `path`, `version`, `cache_dir`, optional
`index_scope`/`provenance`. Each `ld.gwas`/`ld.eqtl` entry names `npz`, `values_key`,
`snp_ids_key`, `positions_key`, `kind`, `panel` and provenance. `ld.heatmap`
explicitly selects `gwas`, `eqtl` or null. Paths resolve relative to the bundle.
No Python code is evaluated from the bundle and no online requests occur.

The composite heatmap first intersects the full reference matrix with the
analysed shared SNP IDs (or actual shared association IDs for no-CS), preserving
original full-matrix indices. Extra reference-only variants are listed as
`omitted_not_shared`. Explicit LD provenance build/chromosome must agree with the
region. Annotation models reject noninteger/nonfinite/reversed coordinates,
invalid strands and child chromosome/parent-bound conflicts. Recognized GENCODE
header releases must agree with the supplied version. A source modified during
index creation is rejected before commit. Supported GFF3 transcript parents are
transcript, mRNA, ncRNA, lnc_RNA, tRNA, rRNA, miRNA, snRNA, snoRNA, scRNA,
primary_transcript, pseudogenic_transcript and processed_transcript; other exon
parent types fail explicitly. A canonical transcript outside the window is named
as unavailable there, without drawing a reversed interval.

The default heatmap uses a restrained continuous white-to-red scale, independently
of the discrete association colors. `heatmap_cmap` (composite) or `cmap`
(standalone helper) accepts another matplotlib colormap. The heatmap **and its
colorbar** remain vector. Unresolved text outside the canvas raises a layout
error before any output is saved; increase width/font spacing or shorten labels.

For an executable, explicitly analytical illustration without a local annotation file, run `python examples/locus_analytical.py --output-prefix figures/analytical`. It creates SVG/PDF/PNG and a metadata sidecar from hand-written inputs, not a biological locus or simulated study.
