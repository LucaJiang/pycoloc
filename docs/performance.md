# Runtime and peak process RSS

## Measured versions and inputs

These existing full-call measurements use pycoloc 0.2.2rc3 and prusie 0.2.3rc4, R 4.4.0 with OpenBLAS 0.3.20, coloc 6.0.3 and susieR 0.16.6. The inference code is unchanged in pycoloc 0.2.2; these are measurements of 0.2.2rc3, not new timings. Each backend uses one numerical thread. The fixed real panel has 100 regions, 200 traits and 103–4162 SNPs per region. It consists of previously selected association windows from one disease/cell-type pairing; it is not a random sample of all loci.

| Task | Comparable cases | Median time, R → Python (ms) | Paired speedup¹ | Faster / tied / slower |
| --- | ---: | ---: | ---: | ---: |
| ABF | 100/100 | 4.33 → 4.29 | 0.928× | 35 / 0 / 65 |
| All component BFs | 99/100 | 40.3 → 13.5 | 2.91× | 99 / 0 / 0 |
| Existing-fit CS pairing | 65/100 | 8.49 → 4.2 | 1.81× | 64 / 0 / 1 |

¹ Geometric mean of per-case median R / median Python runtime. The two displayed time columns are separate medians across cases, so their ratio need not equal the paired summary.

| Task | Median peak process RSS, R → Python (MiB) | Lower / tied / higher in Python |
| --- | ---: | ---: |
| ABF | 259 → 97.7 | 100 / 0 / 0 |
| All component BFs | 260 → 97.5 | 99 / 0 / 0 |
| Existing-fit CS pairing | 260 → 96.9 | 65 / 0 / 0 |

The primary subset requires both calls to succeed and converge, with an available pairing result. Every excluded state remains in the complete table. ABF includes all successful comparisons. Fitting is excluded on both sides of the BF/CS pairing comparisons.

| Task / subset | Count | Speedup median [Q1, Q3]; range | Relative repeat IQR, R / Python (median) |
| --- | ---: | --- | ---: |
| ABF / primary comparable | 100 | 0.939 [0.814, 1.07]; 0.541–1.61 | 0.152 / 0.0579 |
| ABF / all measured | 100 | 0.939 [0.814, 1.07]; 0.541–1.61 | 0.152 / 0.0579 |
| All component BFs / primary comparable | 99 | 2.89 [2.59, 3.23]; 2.11–4.58 | 0.12 / 0.0119 |
| All component BFs / all measured | 100 | 2.88 [2.59, 3.23]; 2.11–4.58 | 0.119 / 0.0119 |
| Existing-fit CS pairing / primary comparable | 65 | 1.89 [1.64, 2.16]; 0.726–2.8 | 0.117 / 0.0465 |
| Existing-fit CS pairing / all measured | 100 | 1.62 [0.217, 2.02]; 0.0874–2.8 | 0.122 / 0.0326 |

Relative repeat IQR is (Q3−Q1)/median across the seven retained times. A runtime ratio below one favors R; ties use exact median equality. These are observed counts without a significance or minimum-percentage threshold.

### Public 500-SNP fixture

| Task | Comparable cases | Median time, R → Python (ms) | Paired speedup¹ | Faster / tied / slower |
| --- | ---: | ---: | ---: | ---: |
| ABF | 2/2 | 4.17 → 4.84 | 0.855× | 0 / 0 / 2 |
| All component BFs | 2/2 | 180 → 48.2 | 3.73× | 2 / 0 / 0 |
| Existing-fit CS pairing | 2/2 | 8.59 → 4.98 | 1.73× | 2 / 0 / 0 |

| Task | Median peak process RSS, R → Python (MiB) | Lower / tied / higher in Python |
| --- | ---: | ---: |
| ABF | 259 → 99.2 | 2 / 0 / 0 |
| All component BFs | 263 → 99 | 2 / 0 / 0 |
| Existing-fit CS pairing | 260 → 97 | 2 / 0 / 0 |

These are the full existing synthetic toy datasets, separately reported from the real panel. D1–D4 are the fine-mapping timing panel; D1/D2 and D3/D4 are the pairing panel. The three additional diagnostic fits are numerically checked but not timed.

### Complete records and counterexamples

The per-case table ([TSV](assets/r_comparison_cases.tsv), [JSON](assets/r_comparison_cases.json)), [all retained repeats](assets/r_comparison_repeats.json), [summary](assets/r_comparison_summary.json) and [parameters/protocol](assets/r_comparison_protocol.json) preserve all workloads and statuses. Per-repeat rows link to separate numerical validation and separate case/backend memory measurements; those quantities are not measured inside every timing repetition. Private genotypes, association files, variant IDs and input paths are not redistributed.

- **ABF:** 65/100 slower primary cases (first IDs: region02, region03, region04, region05, region07, region08, region09, region11); 0/100 higher-RSS cases (first IDs: none). Lowest speed ratio: region49, 0.541×; lowest memory ratio: region24, 2.61×. Paired memory ratio median 2.66 [Q1 2.64, Q3 2.67], range 2.61–2.7; median per-case RSS reduction 62.4%. Summed case medians: R 0.447 s, Python 0.481 s.
- **All component BFs:** 0/99 slower primary cases (first IDs: none); 0/99 higher-RSS cases (first IDs: none). Lowest speed ratio: region08, 2.11×; lowest memory ratio: region48, 2.62×. Paired memory ratio median 2.67 [Q1 2.66, Q3 2.67], range 2.62–2.75; median per-case RSS reduction 62.5%. Summed case medians: R 4.1 s, Python 1.41 s.
- **Existing-fit CS pairing:** 1/65 slower primary cases (first IDs: region45); 0/65 higher-RSS cases (first IDs: none). Lowest speed ratio: region45, 0.726×; lowest memory ratio: region08, 2.58×. Paired memory ratio median 2.68 [Q1 2.67, Q3 2.69], range 2.58–2.71; median per-case RSS reduction 62.7%. Summed case medians: R 0.626 s, Python 0.364 s.

An OS peak includes imports, inputs, the full return and native allocations; it is not incremental algorithm memory. BF/CS consumers load matched complete BF matrices, variant IDs, original CS/member IDs and convergence metadata, without the unused coefficient arrays of an upstream fitting workspace. Fine-mapping retains its full standard result. Runtime and memory conclusions apply to this panel and environment.

All formal calls ran sequentially in one reserved measurement window. The audit recorded 0 pressure-contaminated attempts and 0 process-failed attempts. The predeclared rule retains contaminated attempts and repeats both backends; all original records remain available. Shared-host load monitoring is best effort. The default NumPy coloc path was measured; no fresh Numba compilation or warm-JIT performance claim is made.


## Reproduce a full-call comparison

The [bounded comparison scripts](https://github.com/LucaJiang/pycoloc/tree/main/benchmarks)
run the full `coloc_abf`, `coloc_bf_bf` or `coloc_susie` API against coloc 6.0.3.
Inputs are loaded before the clock; checks, preprocessing and ordinary result
construction remain inside it. Each backend has one warmup and at least seven
retained repeats. Small calls are batched identically. Backend order is alternated
between cases; numerical library threads are fixed to one.

A separate fresh process performs one complete call per case/backend for
**peak process RSS** in MiB. The OS high-water mark includes the interpreter,
imports, loaded inputs, native allocations and the returned object. Numerical
serialization happens after capture. This metric is distinct from SuSiE-RSS,
the statistical summary-data model.

Report median R time / median Python time for each case and paired distributions,
including slower cases, variability and convergence states. Memory ratios are
reported separately; runtime gains do not imply smaller process RSS. Only
measurements with a stated R version, optimized BLAS, thread count and matched
inputs support a comparison with R.

## Implementation scope

The stable H3 numerical reduction uses O(m) work and no m×m configuration matrix.
The full API also validates, aligns and sorts SNPs and constructs pandas tables;
its complexity is not described by that reduction alone. Optional Numba compiles
only the H3 reduction, not input preparation or full result construction. Its
first call includes compilation and must be separated from warm calls.

LD is unnecessary for ABF in either Python or R. Dense LD for upstream
fine-mapping costs O(m²) storage. Plot-only LD thinning leaves analysis inputs
and association tracks unchanged. See [locus plots](locus_plot.md) for rendering
and annotation details, and [agreement](r_agreement.md) for numerical scope.
