# pycoloc

Python Bayesian colocalisation using NumPy/SciPy, with optional Numba for the
linear-time hypothesis kernel. This package implements the statistical scope of
`coloc.abf`, `coloc.bf_bf` and consumption of existing SuSiE fits by `coloc.susie`.
It has no R runtime dependency; ABF and fit consumption do not import prusie.

Author and maintainer: Wenxin Jiang. Release candidate **0.2.2rc3**.

Read the [static documentation](docs/index.html), [500-SNP example guide](docs/example.md), or [executed R agreement report](docs/r_agreement.md). The related Pages site is prebuilt in `docs/`; see the Pages setup instructions in the documentation.

```sh
python -m pip install .
python examples/check_example.py --output-dir example-output
```

The check runs offline once dependencies are installed, using the complete frozen official **synthetic teaching fixture**, without R or prusie. It validates inputs/references by checksum and compares real API calls to separate official-R, mathematical, and native checkpoints. See the guide for the optional native package-pair mode.

Install from this directory with `python -m pip install .` or include optional
JIT and tests with `python -m pip install -c requirements-dev.lock '.[jit,test]'`. Python 3.10+ is required.

```python
from pycoloc import coloc_abf, plot_coloc

# A hand-sized API example, not a biological study or simulated validation.
d1 = dict(snp=['a', 'b', 'c'], position=[100, 200, 300], chromosome='22',
          beta=[0.1, 0.5, -0.2], varbeta=[0.01, 0.02, 0.03], type='cc')
d2 = dict(d1, beta=[0.0, 0.4, 0.1], type='quant', sdY=1.2)
r = coloc_abf(d1, d2, p1=1e-4, p2=1e-4, p12=1e-5,
              prior_weights1={'a': 1, 'b': 2, 'c': 1})
print(r.summary)
print(r.results[['snp', 'SNP.PP.H4', 'SNP.PP.H4.weighted']])
fig, axes = plot_coloc(r, d1, d2, trait1='Case-control trait', trait2='Quantitative trait')
fig.savefig('region.svg')
```

Use signed beta, **SE²** as varbeta, and log(OR) for case-control beta. Quantitative
traits require known sdY or MAF/N for the official estimator. ABF needs no LD.
P-value input additionally requires MAF, N and case fraction s for case-control
traits. When both paths are provided, beta/varbeta takes precedence. Do not use
an assumed sdY=1 unless that trait scale is established.

```python
from pycoloc import coloc_bf_bf, coloc_susie

pair = coloc_bf_bf([1., 4., 0.], [0., 3., 1.],
                   snps1=['a', 'b', 'c'], snps2=['a', 'b', 'c'])
# Given existing native or exported official fits:
# pair = coloc_susie(fit1, fit2, prior_weights1=weights_by_snp)
# fig, axes = plot_coloc(pair, d1, d2, signal_pair=(0, 2))
```

`summary` is always a pandas DataFrame, with one row for ABF or one per retained
component pair. H0 is neither trait associated; H1/H2 only one associated; H3 two
different causal SNPs; H4 a shared causal SNP. Pair-specific H4 values are **not**
a region-level posterior. `results` holds SNP-level probabilities conditional on
H4, meaningful only when H4 itself has support. `priors` records actual priors.
From version 0.1.1 these main SNP fields include both coloc variant weights;
the `.weighted` columns are identical aliases. `summary.h4_status` distinguishes
zero shared evidence (conditional NaN) from floating PP.H4 underflow (conditional
still defined). Partial-zero weights are supported in all three interfaces.
Explicit effect_priors are SDs and require the beta/varbeta branch; P/MAF/N
uses fixed prior scales and rejects explicit overrides.

The source reference is official coloc 6.0.3 commit
`8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45`, selected because the locked CRAN
5.2.3 release lacks variant-specific prior interfaces. This is a focused
implementation, not full R API compatibility. Read
[methods](docs/methods.md), [compatibility](docs/compatibility.md),
[source mapping](docs/reference_mapping.md),
[upstream differences](docs/upstream_differences.md), and
[validation](docs/testing.md) before interpreting weighted comparisons.

GPL-3.0-or-later. See LICENSE, THIRD_PARTY_NOTICES and CITATION.cff.

Publication-size locus figures are available in 0.2.0 via `plot_coloc_locus`,
`prepare_region_annotation`, `plot_gene_track` and `plot_ld_heatmap`. They combine
all H0–H4, two marginal association tracks, local GENCODE genes and real pairwise
LD, with explicit missing-LD colors and vector SVG/PDF plus PNG/JSON export.
Install `python -m pip install -c requirements-dev.lock '.[annotation]'` from this checkout only when indexing local GTF/GFF3. Prepared annotation
models plot offline without that extra. See [the API and local example guide](docs/locus_plot.md).

See [installation](docs/install.md), [all public APIs](docs/api.md), [inputs and priors](docs/inputs.md), [result schema](docs/result_schema.md), [release notes](docs/release_notes.md), and [performance limits](docs/performance.md). The separate repository named **prusie** supplies distribution/import **prusie**. No package-index publication is implied.

See the [name migration and license status](docs/migration.md) for the public rename, saved-result compatibility and retained licensing constraints.

Source repository: [pycoloc](https://github.com/LucaJiang/pycoloc). For Pages setup, see [the guide](docs/pages.md).
Related package: [prusie](https://github.com/LucaJiang/prusie).
