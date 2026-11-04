# Agreement with R

pycoloc 0.2.2 was compared with coloc 6.0.3 on 100 real regions and on the
500-variant D1–D4 toy example from coloc's `coloc_test_data`. Each region has a
GWAS trait and an eQTL trait, with 103–4162 variants. For the fitting comparison,
prusie 0.3.1 was compared with susieR 0.16.6.

The Python results are new. The R results are the stored outputs of an earlier
run (R 4.4.0, OpenBLAS 0.3.20), reused after checking that their inputs match.

## Same inputs

These comparisons give both packages identical inputs, so they test the
colocalisation arithmetic alone. ABF uses the same summary statistics. Bayes-factor
pairing uses all component rows of the same Bayes-factor matrices. Credible-set
pairing uses the same susieR fits and selects their credible-set components.

Each region contributes its largest absolute difference; the table summarises
those values across regions.

| Comparison | Field | Regions | Median | 95th percentile | Maximum |
| --- | --- | ---: | ---: | ---: | ---: |
| ABF | H0–H4 | 100 | 3.3 × 10⁻¹⁶ | 2.7 × 10⁻¹⁵ | 1.1 × 10⁻¹⁴ |
| ABF | `SNP.PP.H4` | 100 | 6.0 × 10⁻¹⁶ | 5.1 × 10⁻¹⁵ | 1.2 × 10⁻¹⁴ |
| All component Bayes factors | H0–H4 | 100 | 8.6 × 10⁻¹⁶ | 1.1 × 10⁻¹⁴ | 1.0 × 10⁻¹¹ |
| All component Bayes factors | `SNP.PP.H4` | 100 | 4.7 × 10⁻¹⁶ | 6.7 × 10⁻¹⁶ | 4.6 × 10⁻¹⁵ |
| Same susieR fits | H0–H4 | 66 | 4.4 × 10⁻¹⁶ | 7.4 × 10⁻¹⁵ | 1.6 × 10⁻¹⁴ |
| Same susieR fits | `SNP.PP.H4` | 66 | 3.5 × 10⁻¹⁶ | 4.8 × 10⁻¹⁶ | 1.4 × 10⁻¹⁵ |

The other 34 regions have no credible set for at least one trait, so there is no
pair to compare; both packages report this. The largest Bayes-factor pairing
difference is in region25 (1.0 × 10⁻¹¹). The largest ABF log Bayes factor
differences are 3.1 × 10⁻¹³ (region01, trait 1) and 5.1 × 10⁻¹³ (region08, trait 2).

## Fitting, then colocalisation

Each trait was fitted with the recorded settings of its case, including the
number of components, sample size and a limit of 100 iterations. The prusie
results went directly into `pycoloc.coloc_susie`; the susieR results went into
`coloc.susie`. Signal pairs are matched by their original component IDs, so the
order in which credible sets are returned does not matter.

Both fits converged in 99 regions. Of these, 65 give at least one pair of
credible sets and 34 have no credible set for at least one trait.

| Comparison | Field | Cases | Median | 95th percentile | Maximum |
| --- | --- | ---: | ---: | ---: | ---: |
| Colocalisation | H0–H4 | 65 regions | 1.8 × 10⁻¹¹ | 4.0 × 10⁻⁹ | 9.7 × 10⁻⁹ |
| Colocalisation | `SNP.PP.H4` | 65 regions | 3.2 × 10⁻¹⁰ | 1.0 × 10⁻⁸ | 7.5 × 10⁻⁸ |
| Fit | `alpha` | 199 traits | 3.8 × 10⁻¹⁰ | 3.7 × 10⁻⁹ | 1.3 × 10⁻⁸ |

The largest differences are in region32 (H0–H4), region08 (`SNP.PP.H4`) and the
GWAS trait of region16 (`alpha`). The two implementations return the same
credible-set components and members in every fitted case, so the remaining
differences come from the separately computed Bayes factors.

The GWAS fit of region25 reaches the 100-iteration limit in both implementations
and is left out of the table above; its largest `alpha` difference is 4.4 × 10⁻¹⁰.

## Toy example

D1/D2 and D3/D4 keep all 500 variants in their original order, and all four fits
converge. Across the two pairs, the largest H0–H4 difference is 9.1 × 10⁻¹⁰ and
the largest `SNP.PP.H4` difference is 3.7 × 10⁻¹⁰. Two extra checks behave the
same in both implementations: a D3 fit stopped after one iteration stays
nonconverged, and a D3 fit with a strict purity threshold returns no credible set,
so its posterior is unavailable rather than zero.

The toy data are simulated, with N = 1000 and stored sdY. Positions are indices,
and genome build, ancestry and effect alleles are not defined. D3 contains two
signals, so ABF on D3 shows software behaviour rather than a valid
single-causal-variant analysis.

## Deliberate differences

For a few inputs pycoloc corrects coloc on purpose. For example, two equal Bayes
factors with prior weights [9, 1] and [1, 1] give conditional shared probabilities
[0.9, 0.1] in pycoloc, while coloc 6.0.3 reports [0.5, 0.5] for its main
`SNP.PP.H4`. These cases are checked against exact calculations instead of R; see
[differences from coloc](upstream_differences.md).

## Reproduce

From a clone, after [installation](install.md):

```sh
python examples/check_example.py --output-dir example-output
pip install prusie==0.3.1
python examples/check_example.py --native-pair --output-dir native-example-output
```

The first command needs neither R nor prusie. The [example guide](example.md)
lists the parameters and where the stored R results come from. Per-field
distributions, every `alpha` layer and the nonconverged cases are in
[numerical_v022.json](assets/numerical_v022.json); the real-region inputs are not
redistributed. Runtime and memory are on the [performance](performance.md) page.
