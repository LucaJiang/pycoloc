# API reference

Annotation serialization helpers and model containers are in `pycoloc.annotation`; inference and plotting entry points are top-level exports.

The following signatures and parameter descriptions describe the public API. See the input and result guides for alignment and interpretation.

## coloc_abf

```python
coloc_abf(dataset1: Mapping, dataset2: Mapping, *, p1: float=0.0001, p2: float=0.0001, p12: float=1e-05, prior_weights1=None, prior_weights2=None, MAF=None, effect_priors: Mapping | None=None, engine: str='numpy') -> ColocResult
```

Colocalise two traits under at most one causal variant per trait.

### Parameters
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

### Returns
ColocResult
    One-row summary, per-SNP ABFs and conditional shared probabilities.

### Raises
ValueError
    Invalid statistics, priors or weights. No overlap instead returns
    status='no_overlap' with an empty summary, never H4=0.

### Notes
SNP.PP.H4 is conditional on H4 and includes both coloc weights;
SNP.PP.H4.weighted is an alias. SNPs are aligned by ID; conflicting
shared positions raise rather than silently changing the analysis support.

### References
coloc 6.0.3, claudia.R: process.dataset, combine_abf_weighted, coloc.abf.

### Examples


```pycon
>>> d = dict(snp=['a', 'b'], beta=[0.2, 0.0], varbeta=[0.01, 0.01], type='cc')
>>> coloc_abf(d, d).nsnps
2
```

## coloc_bf_bf

```python
coloc_bf_bf(bf1, bf2, *, snps1, snps2, p1: float=0.0001, p2: float=0.0001, p12: float=5e-06, prior_weights1=None, prior_weights2=None, null_index1: int | None=None, null_index2: int | None=None, overlap_min: float=0.5, trim_by_posterior: bool=True, engine: str='numpy') -> ColocResult
```

Colocalise every pair of natural-log Bayes-factor rows.

### Parameters
bf1, bf2 : array-like, shape (L, m) or (m,)
    Full component rows of natural log BF; -inf denotes zero BF. An
    explicit final null column changes shape to (L, m+1). No PIPs.
snps1, snps2 : sequence of str
    Unique biological SNP IDs in each matrix's column order, excluding null.
p1, p2, p12 : float
    Positive scalar priors, defaults 1e-4, 1e-4, 5e-6. BF hypothesis priors
    are not adjusted for region size by the upstream function.
prior_weights1, prior_weights2 : mapping, Series or array, optional
    Nonnegative finite biological-SNP weights, bound by ID (arrays bind
    to supplied SNP order). Each side needs positive total on common SNPs.
    Unlike R's matrix-column length check, null weights are not requested.
null_index1, null_index2 : int, optional
    Must equal m for an explicit final null. Null is never inferred.
overlap_min : float
    Minimum fraction of each component's full biological posterior support
    represented by common SNPs. Coloc weights do not change this check.
trim_by_posterior : bool
    Drop insufficient pairs. If all fail, retain reference-compatible rows
    with H1–H4 missing and H0=min(null posterior1, null posterior2).
engine : {'numpy', 'numba'}
    O(m) stable H3 calculation without a SNP-by-SNP configuration matrix.

### Returns
ColocResult
    First trait row varies fastest, matching R expand.grid. idx1/idx2 are
    zero based, idx1_r/idx2_r are R indices. SNP.PP.H4.abf (one pair) or
    SNP.PP.H4.rowK (multiple) gives P(shared SNP | H4, data), including
    both coloc weights. Parallel .weighted fields are identical aliases.

### Raises
ValueError
    Invalid BF, identifiers, priors, weights or null convention.

### Notes
Common columns and weights are always matched by SNP ID. No-overlap,
insufficient-overlap and partial trimming are distinct result statuses.

### References
coloc 6.0.3 susie.R: coloc.bf_bf and logbf_to_pp.

### Examples


```pycon
>>> r = coloc_bf_bf([1., 3.], [2., 4.], snps1=['a', 'b'], snps2=['a', 'b'])
>>> r.summary.idx1.tolist()
[0]
```

## coloc_susie

```python
coloc_susie(fit1, fit2, *, p1: float=0.0001, p2: float=0.0001, p12: float=5e-06, prior_weights1=None, prior_weights2=None, overlap_min: float=0.5, trim_by_posterior: bool=True, engine: str='numpy') -> ColocResult
```

Consume existing SuSiE fits and pair rows selected by actual cs_index.

### Parameters
fit1, fit2 : mapping or object with as_dict()
    Fields: variant_ids (m biological IDs), lbf_variable (L,m or L,m+1),
    sets.cs, sets.cs_index (zero-based original component IDs), converged,
    and optional null_index=m. Every SNP column is retained, including
    SNPs outside credible sets. Fine-mapping priors are not reapplied.
p1, p2, p12 : float
    Positive scalar coloc priors, defaults 1e-4, 1e-4, 5e-6.
prior_weights1, prior_weights2 : mapping, Series or array, optional
    Nonnegative finite coloc weights on biological SNPs of each fit,
    with positive total on each side after intersection.
overlap_min : float
    Posterior overlap threshold, default 0.5.
trim_by_posterior : bool
    Apply full-support overlap trimming, default True.
engine : {'numpy', 'numba'}
    Stable hypothesis kernel implementation.

### Returns
ColocResult
    Pair-specific hypotheses. Original zero-based component IDs and R
    index mapping are explicit. No CS returns status='no_cs'.

### Raises
ValueError
    Invalid result contract or input, including invalid CS component IDs.

### Notes
Unconverged fits produce a warning and metadata, preserving computable
reference behavior. No fitting, R calls, or prusie import occurs here.

### References
coloc 6.0.3 susie.R: coloc.susie and coloc.bf_bf.

### Examples


```pycon
Pass the SusieResult returned by prusie.susie_rss directly as either fit.
```

## plot_coloc

```python
plot_coloc(result: ColocResult, dataset1: Mapping, dataset2: Mapping, *, signal_pair: tuple[int, int] | None=None, trait1: str='Trait 1', trait2: str='Trait 2', figsize: tuple[float, float]=(10, 7), fontsize: float=10, point_size: float=14, annotate_leads: bool=True, highlight_cs: bool=True)
```

Plot two regional marginal-association tracks and exact result summaries.

### Parameters
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

### Returns
figure, axes
    Editable matplotlib figure and two axes; never calls show().

### Raises
ValueError
    Missing coordinates/statistics, chromosome conflicts, unknown signal
    pair, or plotted datasets missing analyzed SNPs.

### Notes
These are marginal tracks, not conditional or residual associations. H0–H4
come directly from a single result row and are never aggregated over pairs.
Save with figure.savefig('region.svg') or figure.savefig('region.png').

### References
The statistical rows use the coloc 6.0.3 conventions documented by pycoloc.

### Examples


```pycon
>>> # fig, axes = plot_coloc(result, d1, d2, signal_pair=(0, 2))
```

## plot_coloc_locus

```python
plot_coloc_locus(result, gwas, eqtl, *, region, annotation=None, target_gene=None, signal_pair=None, gwas_ld=None, eqtl_ld=None, heatmap_ld=None, leads=None, lead_rule='min_p', trait1='GWAS', trait2='eQTL', chromosome_policy='chr', output_prefix=None, formats=('svg', 'pdf', 'png'), figsize=(180 / 25.4, 6.6), fontsize=7.5, point_size=7, max_labels=2, max_heatmap_snps=180, target_mode='union', max_gene_lanes=6, heatmap_cmap='white_red', dpi=180)
```

Create a summary, two association tracks, genes and physical LD triangle.

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

## prepare_region_annotation

```python
prepare_region_annotation(source, *, chromosome, start, end, genome_build, version, cache_dir, source_name='GENCODE', chromosome_policy='chr', index_scope=None, source_provenance=None)
```

Index a local GTF/GFF3 and query complete genes overlapping a region.

No network access is performed. Install ``pycoloc[annotation]`` to index;
prepared RegionAnnotation objects need no gffutils to plot. ``index_scope``
lists all permitted chromosomes for a deliberate chromosome subset, otherwise
the file is declared a complete supplied source. Coordinate-window subsets
are not supported. ``source_provenance`` records an upstream source/hash and
deterministic subsetting rule without replacing the actual indexed-file hash.
A changed source/version/build/indexer/schema gets a separate index identity.
Corrupt or mismatched cache files are rebuilt while holding a POSIX file lock.

## plot_gene_track

```python
plot_gene_track(annotation, *, target_gene=None, ax=None, target_mode='union', max_lanes=6, fontsize=7)
```

Draw all overlapping genes; default exons are the transcript union.

``target_mode='canonical'`` prefers MANE_Select, then Ensembl_canonical,
then a clearly recorded representative fallback. Other genes remain unions.
Labels are bounded and collision-filtered with target priority; gene bodies
are retained even when their labels are suppressed. Full coordinates remain
in ``ax.gene_metadata``; only visible segments are clipped to the window.

## plot_ld_heatmap

```python
plot_ld_heatmap(ld, *, start, end, ax=None, max_snps=180, lead=None, colorbar=True, colorbar_ax=None, allowed_ids=None, cmap='white_red')
```

Draw a vector triangle with physical-coordinate cells and audited thinning.

Only a real pairwise matrix can supply a heatmap. Repeated physical positions
retain one deterministic SNP (lead preferred); rank thinning affects only
this panel. Return the axis with ``ld_metadata`` and an optional colorbar.

## annotation_to_dict

```python
annotation_to_dict(annotation)
```

Return a JSON-serializable complete annotation snapshot.

## annotation_from_dict

```python
annotation_from_dict(value)
```

Load a prepared snapshot without gffutils, files, or network requests.


## Data containers

ColocResult fields are documented in [result schema](result_schema.md). Plotting/annotation containers below preserve supplied source metadata and validate their coordinates/LD geometry.

### LDMatrix

| Field | Type | Default |
|---|---|---|
| values | `object` | `required` |
| snp_ids | `object` | `required` |
| positions | `object` | `required` |
| kind | `str` | `required` |
| panel | `str` | `required` |
| provenance | `dict` | `field(default_factory=dict)` |
| tolerance | `float` | `1e-08` |

### LeadLD

| Field | Type | Default |
|---|---|---|
| values | `object` | `required` |
| snp_ids | `object` | `required` |
| lead_id | `str` | `required` |
| kind | `str` | `required` |
| panel | `str` | `required` |
| provenance | `dict` | `field(default_factory=dict)` |
| tolerance | `float` | `1e-08` |

### ExonModel

| Field | Type | Default |
|---|---|---|
| exon_id | `str` | `required` |
| chromosome | `str` | `required` |
| start | `int` | `required` |
| end | `int` | `required` |
| strand | `str` | `'.'` |

### TranscriptModel

| Field | Type | Default |
|---|---|---|
| transcript_id | `str` | `required` |
| chromosome | `str` | `required` |
| start | `int` | `required` |
| end | `int` | `required` |
| strand | `str` | `required` |
| exons | `tuple[ExonModel, ...]` | `()` |
| tags | `tuple[str, ...]` | `()` |
| name | `str \| None` | `None` |
| biotype | `str \| None` | `None` |

### GeneModel

| Field | Type | Default |
|---|---|---|
| gene_id | `str` | `required` |
| name | `str` | `required` |
| chromosome | `str` | `required` |
| start | `int` | `required` |
| end | `int` | `required` |
| strand | `str` | `required` |
| transcripts | `tuple[TranscriptModel, ...]` | `()` |
| exons | `tuple[ExonModel, ...]` | `()` |
| biotype | `str \| None` | `None` |

### RegionAnnotation

| Field | Type | Default |
|---|---|---|
| chromosome | `str` | `required` |
| start | `int` | `required` |
| end | `int` | `required` |
| genome_build | `str` | `required` |
| genes | `tuple[GeneModel, ...]` | `()` |
| provenance | `dict[str, Any]` | `field(default_factory=dict)` |
| chromosome_policy | `str` | `'chr'` |
