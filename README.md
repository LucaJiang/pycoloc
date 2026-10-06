# pycoloc

Bayesian colocalisation in Python. Compare two traits using approximate Bayes
factors (ABF), component log Bayes factors, or existing SuSiE fits, and plot the
results. The numerical implementation uses NumPy/SciPy with an optional Numba
kernel. R is not needed at runtime.

## Install and run

```sh
git clone https://github.com/LucaJiang/pycoloc.git
cd pycoloc
python -m venv .venv
. .venv/bin/activate
python -m pip install .
python examples/check_example.py --output-dir example-output
```

Python ≥3.10 is required. See [installation](docs/install.md) for dependencies
and tested platforms. Source installation is the supported download route;
no PyPI release or public wheel download is provided.

The check runs offline after installation on a frozen, explicitly **synthetic
500-SNP teaching fixture**. It compares public API calls with official R outputs
and separate mathematical expectations, prints PASS/FAIL and writes a JSON report.
ABF and pairing existing fits do not require prusie; [native fine-mapping](docs/example.md#optional-native-fitting-with-prusie)
is an optional step.

## Compare two traits

```python
from pycoloc import coloc_abf

# Small analytical illustration; these are not biological data.
d1 = dict(snp=['a', 'b', 'c'], beta=[0.1, 0.5, -0.2],
          varbeta=[0.01, 0.02, 0.03], type='cc')
d2 = dict(d1, beta=[0.0, 0.4, 0.1], type='quant', sdY=1.2)
result = coloc_abf(d1, d2, prior_weights1={'a': 1, 'b': 2, 'c': 1})
print(result.summary)
print(result.results[['snp', 'SNP.PP.H4']])
```

Supply signed beta (log OR for case-control) and **SE²** as `varbeta`.
Quantitative traits require known `sdY` or MAF/N for its estimator; do not assume
unit phenotype variance. P-value inputs require MAF, N and, for case-control
traits, case fraction `s`. ABF needs no LD. Read [inputs and priors](docs/inputs.md)
for alignment and model requirements.

H0 means neither trait is associated; H1/H2 mean only one is associated; H3 means
different causal variants; H4 means a shared causal variant. `SNP.PP.H4` is
conditional on H4, so interpret it alongside H4 support. ABF assumes at most one
causal variant per trait. `coloc_susie` pairs the original credible-set components
of supplied fits; pair-specific posteriors do not define a single region-wide H4.

## Documentation and validation

Read the [documentation site](https://lucajiang.github.io/pycoloc/), or the
Markdown [example](docs/example.md), [API](docs/api.md),
[result schema](docs/result_schema.md) and [locus plotting guide](docs/locus_plot.md).
The [executed R agreement report](docs/r_agreement.md) separates agreement with
coloc 6.0.3 from intentional corrections to weighted SNP probabilities, alignment
and numerical edge cases. It retains the raw R differences and independent
mathematical checks. See [methods](docs/methods.md) and
[limitations](docs/compatibility.md) before interpreting an analysis.

Maintained by Wenxin Jiang. See [contributing](docs/contributing.md),
[changelog](CHANGELOG.md) and [citation](CITATION.cff).
GPL-3.0-or-later; [license](LICENSE) and [third-party notices](THIRD_PARTY_NOTICES)
retain the upstream code and teaching-data terms.
