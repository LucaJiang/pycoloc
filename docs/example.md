# Run the 500-SNP reference check

The example executes `coloc_abf`, `coloc_bf_bf`, and `coloc_susie` on a fixed, redistributable **synthetic upstream teaching fixture**. Each of its four datasets contains 500 distinct SNP-like variables. It is software validation, not evidence from a new biological study.

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

The command prints `PASS` or `FAIL`, field-level maximum absolute errors, and the thresholds used. It writes `report.json` and `actual_results.json` to your output directory and exits with a nonzero code for a substantive failure. A package version change is recorded transparently; its string alone cannot cause a numerical failure. Output inside the frozen data directory is rejected. The runner never updates its expected results.

The [executed R agreement report](r_agreement.md) explains which differences are expected and why. A final `PASS` means the stated checks passed; it does **not** mean every upstream field was identical.

## The data and model

All stored D1, D2, D3 and D4 data were exported from coloc's pinned `coloc_test_data.rda`, unchanged and without outcome-based selection. Each has N=1000, quantitative beta/varbeta and its actual stored sdY. The stored causal labels identify s105 in D1/D2/D4, and s105 plus s89 in D3; these are synthetic variable identities. Positions are indices, not genomic coordinates. Build, ancestry, effect alleles and counting alleles are unavailable and explicitly null. Do not convert these into biological claims.

The bundle supplies numeric float64 NPZ with `allow_pickle=False`, portable TSV/JSON, and a variant manifest that maps both LD axes. LD values are signed correlation r. The full input bundle is about 9 MB, with four complete 500×500 matrices. No human study data or annotation database is included.

The literal default ABF calls pass only the two datasets: p1=p2=1e-4 and p12=1e-5. A second pair of ABF calls explicitly uses p12=5e-6, matching the SuSiE pairing priors. Both use the standard quantitative effect-prior SD `0.15 × sdY`. It requires no LD. Frozen SuSiE fits use known-n z-RSS, `z=beta/sqrt(varbeta)`, L=10, max_iter=100, tol=0.001, coverage=0.95, min_abs_corr=0.5 and full 500-variable purity. Prior variance uses `optim`; residual variance remains fixed. All arguments are recorded in `examples/data/parameters.json`.

Two explicitly labelled diagnostics change one setting: D3_strict_purity uses min_abs_corr=1 and actually returns no CS; D3_one_iteration uses max_iter=1 and actually fails to converge. They do not change the primary fit protocol.

## Optional native fitting with prusie

Install [prusie](https://github.com/LucaJiang/prusie) into the same activated
environment, with its Rust toolchain, C linker and Python headers available.
Starting in the pycoloc checkout:

```sh
git clone https://github.com/LucaJiang/prusie.git ../prusie
python -m pip install ../prusie
```

Then run:

```sh
python examples/check_example.py --native-pair --output-dir native-example-output
```

This additionally calls `prusie.susie_rss` for D1–D4 with the same known-n parameters and passes those native result objects to `coloc_susie`. It checks PIP absolute error ≤1e-5 with zero relative tolerance, input/model and probability validity. It records convergence, iterations and downstream colocalisation differences separately. An end-to-end result combines fine-mapping and colocalisation numerical behavior; it is a distinct comparison from pairing identical frozen fits.

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
