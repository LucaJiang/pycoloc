# pycoloc: Bayesian colocalisation in Python

[![tests](https://github.com/LucaJiang/pycoloc/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/LucaJiang/pycoloc/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)
[![license: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue)](LICENSE)

pycoloc asks whether two traits share a causal variant in a region. It implements
the methods of the R package [coloc](https://github.com/chr1swallace/coloc):
single-variant colocalisation with approximate Bayes factors (ABF), and
colocalisation of SuSiE fine-mapping results, which compares every pair of
credible-set signals. Given the same inputs, it reproduces coloc 6.0.3 to about
10⁻¹⁴.

## Features

- **ABF colocalisation** (`coloc_abf`) from effect estimates and their variances,
  or from P values, MAF and sample size; quantitative and case-control traits; the
  standard p1, p2 and p12 priors, with optional variant-specific prior weights.
- **SuSiE colocalisation** (`coloc_susie`, `coloc_bf_bf`) from the per-component
  log Bayes factors of an existing fit. It accepts [prusie](https://github.com/LucaJiang/prusie)
  results directly, and exported susieR fits through a documented result format.
- **Results** as one `ColocResult`: posterior support for H0–H4 for each signal pair,
  per-variant `SNP.PP.H4`, original credible-set component IDs and convergence
  metadata. Missing credible sets are reported as unavailable rather than H4 = 0.
- **Locus plots** of both association tracks with LD colouring and local gene annotation.
- **Corrections** to a few upstream edge cases, such as variant-weighted `SNP.PP.H4`;
  see [differences from coloc](docs/upstream_differences.md).

## Agreement with R

The comparison uses 100 real regions (103–4162 variants, a GWAS trait and an eQTL
trait in each) and coloc's 500-variant toy example. Each row reports the largest
absolute difference from the R result across all regions.

| Comparison with R | Regions | H0–H4 | `SNP.PP.H4` |
| --- | ---: | ---: | ---: |
| ABF, same summary statistics | 100 | 1.1 × 10⁻¹⁴ | 1.2 × 10⁻¹⁴ |
| Pairing the same component Bayes factors | 100 | 1.0 × 10⁻¹¹ | 4.6 × 10⁻¹⁵ |
| Pairing the same susieR fits | 66 | 1.6 × 10⁻¹⁴ | 1.4 × 10⁻¹⁵ |
| prusie 0.3.1 fit → pycoloc vs susieR 0.16.6 fit → coloc | 65 | 9.7 × 10⁻⁹ | 7.5 × 10⁻⁸ |

The last row includes fitting. Both implementations found the same credible sets
in every region; 34 regions have no credible set for at least one trait, and one
region (the GWAS fit of region25) reaches the 100-iteration limit in both
implementations and is
reported separately. [Agreement with R](docs/r_agreement.md) gives the full distributions.

## Performance

Runtime and peak memory against coloc 6.0.3 in R 4.4.0 with OpenBLAS, one thread,
on the same 100 regions. These were measured with pycoloc 0.2.2rc3; the inference
code has not changed since.

| Task | Regions | Median time, R → pycoloc | Speedup¹ | Regions faster | Median peak RSS, R → pycoloc |
| --- | ---: | ---: | ---: | ---: | ---: |
| ABF | 100 | 4.33 → 4.29 ms | 0.93× | 35 | 259 → 98 MiB |
| Pairing all component Bayes factors | 99 | 40.3 → 13.5 ms | 2.91× | 99 | 260 → 98 MiB |
| Pairing credible sets of existing fits | 65 | 8.49 → 4.20 ms | 1.81× | 64 | 260 → 97 MiB |

¹ Geometric mean of the per-region R/pycoloc time ratio.

ABF takes a few milliseconds in both and is slightly slower in pycoloc for most
regions. Peak RSS is measured for a fresh process, so it mostly reflects the
Python and R runtimes rather than the algorithm. See [performance](docs/performance.md).

## Installation

pycoloc requires Python 3.10 or newer. Install the 0.2.2 wheel from the GitHub Release:

```sh
pip install https://github.com/LucaJiang/pycoloc/releases/download/v0.2.2/pycoloc-0.2.2-py3-none-any.whl
```

or install from the tagged source:

```sh
pip install "git+https://github.com/LucaJiang/pycoloc.git@v0.2.2"
```

pycoloc is not on PyPI. The PyPI package named `pycoloc` is a different project,
so `pip install pycoloc` will not install this one. To fit SuSiE models in Python,
also run `pip install prusie`. See [installation](docs/install.md) for the optional
Numba and annotation extras.

## Quick start

The example data are coloc's simulated `coloc_test_data` (four quantitative traits,
500 variants), included in the repository. Clone it and run from the checkout:

```python
import json
from pathlib import Path

import numpy as np
from pycoloc import coloc_abf

data = Path("examples/data")
metadata = json.loads((data / "metadata.json").read_text())
with np.load(data / "inputs.npz", allow_pickle=False) as arrays:
    datasets = [
        dict(
            snp=metadata[key]["snp"],
            type=metadata[key]["type"],
            N=metadata[key]["N"],
            sdY=metadata[key]["sdY"],
            beta=arrays[key + "_beta"],
            varbeta=arrays[key + "_varbeta"],
        )
        for key in ("D1", "D2")
    ]

result = coloc_abf(*datasets)
print(result.summary[[f"PP.H{i}.abf" for i in range(5)]])
print(result.results[["snp", "SNP.PP.H4"]].nlargest(5, "SNP.PP.H4"))
```

H0 means neither trait has an association in the region, H1 and H2 that only one
does, H3 that both do with different causal variants, and H4 that they share one.
Here `PP.H4.abf` ≈ 0.999 and s105 has the highest `SNP.PP.H4`. `SNP.PP.H4` ranks
variants given that H4 is true; it is not a fine-mapping PIP.

To colocalise SuSiE fits instead, fit each trait once with prusie and pass the
results to `coloc_susie`; the [example](docs/example.md#fit-once-and-reuse-the-result)
shows this on the same data.

## Which function to use

| Input | Function | Assumption |
| --- | --- | --- |
| Effect estimates and variances, or P values with MAF and N | `coloc_abf` | At most one causal variant per trait; no LD needed |
| Per-component log Bayes factors | `coloc_bf_bf` | One posterior per pair of components |
| Existing SuSiE fits | `coloc_susie` | One posterior per pair of credible sets |

With several signals, each row's H4 describes one pair of signals, not the whole
region. SuSiE colocalisation needs signed LD that matches the effect alleles, and
fits that converged; see [inputs](docs/inputs.md).

## Documentation

The [documentation site](https://lucajiang.github.io/pycoloc/) covers the
[API](docs/api.md), [result tables](docs/result_schema.md), [locus plots](docs/locus_plot.md)
and [methods](docs/methods.md). To check an installation against the stored R
results, run `python examples/check_example.py --output-dir example-output`
from the checkout.

## Citation and license

Please cite the original coloc papers listed in [methods](docs/methods.md) and the
pycoloc version you used ([CITATION.cff](CITATION.cff)). pycoloc is maintained by
Wenxin Jiang and licensed under GPL-3.0-or-later ([LICENSE](LICENSE)); upstream
attribution and example-data terms are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
