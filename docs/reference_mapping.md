# Pinned source mapping

Official source: [pinned coloc 6.0.3](https://github.com/chr1swallace/coloc/tree/8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45).

| Python location | Pinned R function | Statistical step |
|---|---|---|
| abf.process_dataset | claudia.R process.dataset | beta/variance priority; cc/quant branches |
| abf.process_dataset | approx.bf.estimates, approx.bf.p, Var.data, Var.data.cc | Wakefield ABF and scales |
| abf.process_dataset | sdY.est | Regression through origin for phenotype variance |
| _validation | check.R check_dataset; claudia.R adjust_prior | Shape/domain checks; region-size adjustment |
| _math.combine | combine.abf, combine_abf_weighted | H0–H4 and Q-scaled weights; corrected main shared conditional from the same q12 |
| _math.offdiag_logsum | logdiff(logsum(l1)+logsum(l2),logsum(l1+l2)) | Algebraically identical H3 via stable positive summation |
| abf.coloc_abf | claudia.R coloc.abf | Full-support priors, intersection, per-SNP result |
| pairing._posterior_support | susie.R logbf_to_pp | Conditional biological support and null posterior |
| pairing.coloc_bf_bf | susie.R coloc.bf_bf | Pair order, overlap trim, hit and summary rows |
| pairing.coloc_susie | susie.R coloc.susie | Actual CS component rows, complete BF support |
| plotting.plot_coloc | Original implementation | Regional marginal tracks with exact result annotations and corrected-main shared candidate |

Correctness checkpoint 0.1.1 / result schema 2: `_math.h4_status` distinguishes
true zero log evidence from probability underflow; summary logU3/logU4 retain
diagnostic evidence. `_validation.weights` admits partial zeros uniformly.
`abf.process_dataset` explicitly rejects effect_priors on P/MAF/N, retaining
the fixed `approx.bf.p` scales. Quantitative beta prior SD is multiplied by sdY,
then squared, as the pinned `approx.bf.estimates` code actually does (despite
its parameter documentation calling the inputs variances).
Counterexamples and regression tests are in [upstream differences](upstream_differences.md).

The R source was read locally from the locked archive, not a floating branch or
third-party Python implementation. The source archive SHA256 is
f6580ab471e1bd9f3acc45e23a73bb2ee16fe3c0a45e01cd8cc008854661b479.
The current toy-example fit exports use [susieR 0.16.6](https://github.com/stephenslab/susieR/tree/8e56a8e038e989856d106d9ca5175cc664fea9d2), commit
`8e56a8e038e989856d106d9ca5175cc664fea9d2`. Its exact archive/data hashes and
R/dependency environment are in `examples/data/reference.lock.json` and
`reference_environment.json`. Methods cited by upstream: Giambartolomei et al., PLoS Genetics,
doi:10.1371/journal.pgen.1004383; Wallace, PLoS Genetics,
doi:10.1371/journal.pgen.1008720 and doi:10.1371/journal.pgen.1009440.
These are upstream method citations, not identifiers for this Python package.

## Locus feature 0.2.0

`annotation.models` contains original provider-independent models, inclusive exon
union and documented transcript selection. `annotation.local` adapts official
gffutils SQLite indexing/overlap/parent queries with local provenance and locking.
`locus_ld` implements explicit LD validation, ID mapping and audited plot-only
thinning. `locus` implements the five-panel layout and physical triangle; these
are new visualization APIs, not ports of statistical R functions. See [locus plots](locus_plot.md) and the annotation/locus tests for these interfaces.
