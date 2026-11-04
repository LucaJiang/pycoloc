# Numerical agreement with R

pycoloc **0.2.2** was evaluated against coloc **6.0.3**, using identical inputs
for 100 frozen real regions and the synthetic 500-SNP D1–D4 toy example from
coloc's `coloc_test_data`. The native fitting comparison uses **prusie 0.3.1**
from PyPI against susieR **0.16.6**.

The Python calculations were rerun. R outputs were reused after verification of
the recorded input hashes; their recorded environment is R 4.4.0 with OpenBLAS
0.3.20. No new R fitting, runtime or memory measurements were made.

## Identical-input comparisons

These comparisons separate colocalisation arithmetic from differences in fitting.
ABF uses the same association inputs. Full-BF pairing uses all supplied component
rows; existing-fit pairing selects original credible-set component IDs from the
same official R fits and retains every SNP column.

The table summarizes the **largest absolute difference within each region**,
then takes the median, 95th percentile and maximum across those regional values.
All 100 regions are accounted for; no-CS results have no numeric posterior rows.

| Comparison / field | Cases with values | Median | 95th percentile | Maximum |
| --- | ---: | ---: | ---: | ---: |
| ABF: H0-H4 | 100 | 3.33e-16 | 2.74e-15 | 1.07e-14 |
| ABF: SNP.PP.H4 | 100 | 5.97e-16 | 5.12e-15 | 1.24e-14 |
| Full component BFs: H0-H4 | 100 | 8.64e-16 | 1.07e-14 | 1.02e-11 |
| Full component BFs: SNP.PP.H4 | 100 | 4.72e-16 | 6.73e-16 | 4.55e-15 |
| Same official fits: H0-H4 | 66 | 4.44e-16 | 7.36e-15 | 1.6e-14 |
| Same official fits: SNP.PP.H4 | 66 | 3.47e-16 | 4.82e-16 | 1.44e-15 |

The largest full-BF hypothesis difference is in region25 (1.02×10⁻¹¹).
The largest ABF logBF differences are 3.13×10⁻¹³ for trait 1 (region01) and
5.12×10⁻¹³ for trait 2 (region08). With the same inputs and component rows,
these differences reflect floating-point evaluation of the same calculations.

## Native fitting followed by colocalisation

Each trait was fitted with the recorded case parameters, including the recorded
component count, sample size and a 100-iteration limit. Native fits were passed
directly to `pycoloc.coloc_susie`; the reference pairs official susieR fits through
`coloc.susie`. Component IDs determine pair identity; CS return order is immaterial.

Among 100 real regions, 99 have mutually converged fits: 65 return posterior
pairs and 34 have no CS for at least one trait. The 199 mutually converged trait
fits supply the alpha comparison below.

| Comparison / field | Cases with values | Median | 95th percentile | Maximum |
| --- | ---: | ---: | ---: | ---: |
| Native fit → coloc: H0-H4 | 65 | 1.77e-11 | 4e-09 | 9.7e-09 |
| Native fit → coloc: SNP.PP.H4 | 65 | 3.23e-10 | 1.03e-08 | 7.47e-08 |
| Native fit: alpha | 199 | 3.79e-10 | 3.68e-09 | 1.28e-08 |

The largest hypothesis difference occurs in region32; the largest conditional
SNP difference occurs in region08. The largest alpha difference is in
region16/gwas. CS component IDs and membership agree for all 206 fitted cases,
including the toy diagnostics. These downstream differences arise when separately
fitted component Bayes factors enter pairing; they do not reflect different CS
membership on this panel.

**region25/gwas** did not converge in either implementation after 100 iterations.
Its region is recorded separately and is excluded from the native distributions
above. The intentional **D3_one_iteration** toy diagnostic also remains
nonconverged in both implementations after one iteration.

## The 500-SNP toy example

D1/D2 and D3/D4 retain all 500 SNPs in their original order. All four primary
native fits converge. Across the two primary native pairs, the largest absolute
H0–H4 difference is 9.14×10⁻¹⁰ and the largest conditional SNP difference is
3.66×10⁻¹⁰. The strict-purity D3 diagnostic returns no CS in both implementations;
its posterior is unavailable, not zero.

The toy data are synthetic, with stored N=1000 and sdY. Positions are indices;
genome build, ancestry and effect alleles are unavailable. D3 contains two signals,
so its ABF comparison illustrates software behavior rather than validating a
single-causal-variant model for those data.

## Intentional differences from upstream

Corrected weighted and partial-zero behavior is checked against independent
mathematical expectations. For example, two equal BFs with weights [9,1] and
[1,1] yield conditional shared probabilities [0.9,0.1]; pinned R's main SNP
posterior is [0.5,0.5]. The full explanation and other scientific corrections
remain in [upstream differences](upstream_differences.md).

The numerical scales reported here describe reproducibility of these computations;
they are not statistical evidence for colocalisation. Interpret H4 alongside the
model assumptions, allele/LD alignment and convergence, and interpret SNP.PP.H4
conditional on H4.

## Reproduce the bundled comparisons

After [installation](install.md), run from the checkout:

```sh
python examples/check_example.py --output-dir example-output
python -m pip install prusie==0.3.1
python examples/check_example.py --native-pair --output-dir native-example-output
```

The first command needs neither R nor prusie. See the [example guide](example.md)
for parameters and frozen reference provenance, and [comparison scripts](https://github.com/LucaJiang/pycoloc/tree/main/benchmarks)
for the full-BF interface. [Numerical distributions and largest cases](assets/numerical_v022.json)
include every field, original alpha layer and individually identified
nonconverged case without redistributing private study inputs.

[Runtime and peak process RSS](performance.md) retain their original measured
pycoloc 0.2.2rc3 identity. Inference code is unchanged in 0.2.2.
