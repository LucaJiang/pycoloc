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
fits above converge on this fixture; inspect convergence and LD consistency for
new data. Both LD axes must match the signed effects' SNP and allele order.
Dense float64 LD stores about 8m² bytes per m-SNP matrix before fitting workspaces.
There is no automatic allele harmonization or LD repair.

Each summary row names original zero-based component IDs. Read its
`snp_pp_column` to find the corresponding conditional SNP posterior. Missing CS
returns `status='no_cs'`; it is unavailable evidence, not H4=0. Existing official
R fits can also be consumed through the [documented result contract](compatibility.md).

## Run offline after installation

Follow [installation](install.md) to clone the repository, create an environment
and install the package. From that checkout:

```sh
python examples/check_example.py --help
python examples/check_example.py --output-dir example-output
```

Installation may download dependencies. Once installed, the check runs offline
without R or external study data. It imports the installed package and uses the
frozen files in `examples/data/`.

## What the check does

1. Verify SHA256 checksums of the input, official-reference and native-checkpoint files.
2. Validate the 500 unique IDs, effect/variance arrays, signed LD order, symmetry and diagonal.
3. Execute unweighted quantitative ABF comparisons for D1/D2 and D3/D4 against official coloc 6.0.3, both with literal public defaults and with the separately recorded common prior settings.
4. Consume the frozen official SuSiE fits and compare component-pair H0–H4, SNP posteriors, IDs, leads and probability validity. Call BF pairing separately on the same full component rows.
5. Check no-CS as unavailable and propagate an explicitly nonconverged frozen fit with its warning.
6. Check corrected weighted and partial-zero behavior against independently enumerated mathematical expectations, preserving the raw R differences as diagnostics.
7. Compare against a separate frozen pycoloc regression checkpoint. It is labelled as native output and never used as an independent R or mathematical oracle.

The command writes field-level numerical differences to `report.json` and results to
`actual_results.json`. It exits nonzero for an unsuccessful software check and never
updates frozen expectations. The [R agreement report](r_agreement.md) explains the
observed differences, including the intentional weighted-field corrections.

## The data and model

All stored D1, D2, D3 and D4 data were exported from coloc's pinned `coloc_test_data.rda`, unchanged and without outcome-based selection. Each has N=1000, quantitative beta/varbeta and its actual stored sdY. The stored causal labels identify s105 in D1/D2/D4, and s105 plus s89 in D3; these are synthetic variable identities. Positions are indices, not genomic coordinates. Build, ancestry, effect alleles and counting alleles are unavailable and explicitly null. Do not convert these into biological claims.

The bundle supplies numeric float64 NPZ with `allow_pickle=False`, portable TSV/JSON, and a variant manifest that maps both LD axes. LD values are signed correlation r. The full input bundle is about 9 MB, with four complete 500×500 matrices. No human study data or annotation database is included.

The literal default ABF calls pass only the two datasets: p1=p2=1e-4 and p12=1e-5. A second pair of ABF calls explicitly uses p12=5e-6, matching the SuSiE pairing priors. Both use the standard quantitative effect-prior SD `0.15 × sdY`. It requires no LD. Frozen SuSiE fits use known-n z-RSS, `z=beta/sqrt(varbeta)`, L=10, max_iter=100, tol=0.001, coverage=0.95, min_abs_corr=0.5 and full 500-variable purity. Prior variance uses `optim`; residual variance remains fixed. All arguments are recorded in `examples/data/parameters.json`.

Two explicitly labelled diagnostics change one setting: D3_strict_purity uses min_abs_corr=1 and actually returns no CS; D3_one_iteration uses max_iter=1 and actually fails to converge. They do not change the primary fit protocol.

## Optional native fitting with prusie

Install prusie 0.3.1 from PyPI into the same environment:

```sh
python -m pip install prusie==0.3.1
```

Then run:

```sh
python examples/check_example.py --native-pair --output-dir native-example-output
```

This additionally calls `prusie.susie_rss` for D1–D4 with the same known-n parameters and passes those native result objects to `coloc_susie`. It compares PIPs with the frozen R fits and checks input, model and probability validity. It records convergence, iterations and downstream colocalisation differences separately. An end-to-end result combines fine-mapping and colocalisation numerical behavior; it is a distinct comparison from pairing identical frozen fits.

The optional import is made only when this flag is present. pycoloc does not install or invoke prusie implicitly. Consult the agreement report for the actual pair of package versions tested.

## Frozen reference provenance and regeneration

The bundle's `README.md`, `LICENSE.GPL-3`, `reference.lock.json`, `reference_environment.json` and `checksums.json` contain attribution, redistribution terms, pinned archives/data hashes, and actual R/dependency versions. coloc is 6.0.3 at commit `8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45`; susieR is 0.16.6 at `8e56a8e038e989856d106d9ca5175cc664fea9d2`. No floating R reference is used.

Reference regeneration is optional and writes a new directory. Follow the exact recipe in `examples/data/README.md` with the locked R 4.4.0 environment. The original generator is preserved for inspection and is **not run**. The separate `examples/data/default_abf/` supplement contains its own frozen R results, checksums and regeneration script for the two literal default-prior calls. Verify the original RDA hash, execute `export_reference.R`, then convert its output with `convert_exports.py`. Never replace the shipped expectations as a way to make a failed check pass.

## Deliberate failure tests

From the checkout, install the test extra and run:

```sh
python -m pip install -c requirements-dev.lock '.[test]'
python -m pytest tests/test_example.py -q
```

Tests run the CLI with no executable search path, reject a changed input checksum, and alter an expected H4 value while repairing that file's checksum to verify that a real numerical mismatch also fails. A separate import/process guard proves the standalone numerical path does not import prusie/rpy2 or spawn an external computation.
