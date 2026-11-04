# Analyze and validate the 500-SNP example

Analyze a complete **synthetic 500-SNP toy example** from official coloc,
then check the results against frozen R outputs. The four bundled datasets
support ABF analysis, existing-fit pairing and optional native fine-mapping.

## Analyze the data

After [installation](install.md), run `python examples/quickstart.py` from the
checkout. The script loads both complete 500-SNP quantitative datasets, calls
`coloc_abf` and prints the hypothesis posteriors and leading conditional SNP
probabilities. The [overview](index.md#quick-start) shows the code and its interpretation.

## Fit once and reuse the result

When several traits or genes share a region, fit each one once at the intended
SNP support and settings. Keep those fits for subsequent pairing. After installing
prusie, this example uses the same stored D1/D2 arrays as the ABF quick start:

```python
import json
from pathlib import Path
import numpy as np
from prusie import susie_rss
from pycoloc import coloc_susie

data = Path("examples/data")
metadata = json.loads((data / "metadata.json").read_text())
fits = {}
with np.load(data / "inputs.npz", allow_pickle=False) as arrays:
    for key in ("D1", "D2"):
        fits[key] = susie_rss(
            z=arrays[key + "_beta"] / np.sqrt(arrays[key + "_varbeta"]),
            R=arrays[key + "_LD"],
            n=metadata[key]["N"],
            variant_ids=metadata[key]["snp"],
            L=10,
            max_iter=100,
            tol=0.001,
            coverage=0.95,
            min_abs_corr=0.5,
            n_purity=500,
            estimate_residual_variance=False,
        )

print({key: fit.converged for key, fit in fits.items()})
paired = coloc_susie(fits["D1"], fits["D2"])
print(paired.summary[["idx1", "idx2", "PP.H3.abf", "PP.H4.abf"]])
```

Reuse `fits` for other pairings without calling `susie_rss` again. The two native
fits above converge on these data; inspect convergence and LD consistency for
new data. Both LD axes must match the signed effects' SNP and allele order.
Dense float64 LD stores about 8m² bytes per m-SNP matrix before fitting workspaces.
There is no automatic allele harmonization or LD repair.

Each summary row names original zero-based component IDs. Read its
`snp_pp_column` to find the corresponding conditional SNP posterior. Missing CS
returns `status='no_cs'`; it is unavailable evidence, not H4=0. Existing official
R fits can also be consumed through the [documented result contract](compatibility.md).

## Check an installation against R

From a clone, after [installation](install.md):

```sh
python examples/check_example.py --output-dir example-output
```

The check runs offline, without R, using the installed package and the stored
files in `examples/data/`. It:

1. verifies the SHA256 checksums of the inputs and stored results;
2. checks the 500 variant IDs, the effect and variance arrays, and that each LD
   matrix is symmetric with a unit diagonal and matches the variant order;
3. runs ABF for D1/D2 and D3/D4 and compares with coloc 6.0.3, once with the
   default priors and once with p12 = 5 × 10⁻⁶ (the SuSiE default);
4. pairs the stored susieR fits with `coloc_susie`, and their full component rows
   with `coloc_bf_bf`, and compares H0–H4, `SNP.PP.H4`, component IDs and leads;
5. checks that a fit without credible sets is reported as unavailable and that a
   nonconverged fit is passed through with a warning;
6. checks the corrected weighted and partial-zero cases against exact calculations,
   keeping the raw R values for comparison;
7. compares with a stored pycoloc result, as a regression check.

It writes the per-field differences to `report.json` and the results to
`actual_results.json`, exits with a nonzero status if a check fails, and never
rewrites the stored results. [Agreement with R](r_agreement.md) discusses the
differences.

## The data

D1–D4 are exported unchanged from coloc's `coloc_test_data.rda`. Each has
N = 1000, quantitative beta/varbeta and its stored sdY. The simulated causal
variant is s105 in D1, D2 and D4, and s105 and s89 in D3. Positions are indices,
not genomic coordinates, and genome build, ancestry and effect alleles are not
defined.

The arrays are stored as float64 NPZ (loaded with `allow_pickle=False`), with
TSV/JSON copies and a variant manifest for both LD axes. LD is the signed
correlation r. The bundle is about 9 MB, mostly the four 500 × 500 LD matrices.

ABF uses the default priors p1 = p2 = 10⁻⁴ and p12 = 10⁻⁵, plus a second run with
p12 = 5 × 10⁻⁶, and the quantitative effect prior SD 0.15 × sdY. The stored SuSiE
fits use summary statistics with known n, `z = beta / sqrt(varbeta)`, L = 10,
max_iter = 100, tol = 0.001, coverage = 0.95, min_abs_corr = 0.5 and purity over
all 500 variants; the prior variance is estimated with `optim` and the residual
variance is fixed. `examples/data/parameters.json` records every argument.

Two extra D3 fits each change one setting: min_abs_corr = 1 gives no credible
set, and max_iter = 1 gives a nonconverged fit.

## Fitting with prusie

```sh
pip install prusie==0.3.1
python examples/check_example.py --native-pair --output-dir native-example-output
```

This also fits D1–D4 with `prusie.susie_rss` using the same settings, compares the
PIPs with the stored susieR fits, and passes the prusie results to `coloc_susie`.
Convergence, iteration counts and the colocalisation differences are reported
separately. prusie is imported only when `--native-pair` is given.

## Where the stored results come from

`examples/data/` contains the attribution and licence (`README.md`,
`LICENSE.GPL-3`), the pinned sources and data hashes (`reference.lock.json`,
`checksums.json`) and the R environment (`reference_environment.json`). The
references are coloc 6.0.3 at commit `8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45`
and susieR 0.16.6 at commit `8e56a8e038e989856d106d9ca5175cc664fea9d2`.

To regenerate the references into a new directory, follow `examples/data/README.md`
with R 4.4.0: verify the hash of the original RDA file, run `export_reference.R`
and convert its output with `convert_exports.py`. `examples/data/default_abf/`
holds the separate results and script for the two default-prior ABF calls.

## Tests of the check itself

```sh
pip install '.[test]'
pytest tests/test_example.py -q
```

These tests change an input file and, separately, an expected H4 value (with a
matching checksum) and confirm that the check fails in both cases. They also
confirm that the standalone check neither imports prusie or rpy2 nor starts an
external process.
