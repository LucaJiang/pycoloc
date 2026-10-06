# 0.2.2rc3 (2026-10-06)

- Update optional fitting integration, runtime labels, tests and docs to import prusie.
- Preserve statistical code, input/reference bytes and existing licenses.
- See [migration notes](docs/migration.md). Earlier entries retain prior-version history.

# 0.2.2rc2

- Add offline public teaching-fixture checks with separate official-R and mathematical expectations.
- Add static Pages documentation, executed agreement report and local integrity/link checks.
- Preserve inference defaults and existing corrected-field semantics.

# Changelog

## 0.2.2rc1

- Preserve the current corrected colocalisation, locus plotting and local annotation implementation.
- Add independent repository packaging, API/schema documentation, examples and public CI.

This local export has no assigned remote, DOI or package-index release. No publication or default-environment replacement was performed.

The development starting point was experimental pyrsusie0.2.2rc17; the separately accepted0.2.0 installation remains unchanged. Python/Cargo/CITATION versions identify this cleaned source.

The validation protocol uses all100 frozen regions/200 trait cases, L=5, max_iter=100 and tol=.001 against pinned susieR0.16.6. PIP absolute error ≤1e-5 (rtol0), probability validity and unchanged input/model semantics define numerical acceptance. Intermediate/credible-set and downstream colocalisation changes remain diagnostics. One known GWAS case reaches100 iterations without convergence and must remain visible in the accompanying validation report.

Local checks use CPython3.12/Linux x86-64. Public tests and examples are self-contained analytical inputs; private benchmark input data and annotation databases are not redistributable artifacts here. The source includes a future GitHub CI workflow; no remote CI success is claimed. Formal measurements and local package/check results accompany the delivery separately. No universal speed or hardware bottleneck claim is made.

## Earlier semantic changes

Corrected SNP.PP.H4 weighting, ID-bound prior alignment, explicit effect-prior SD propagation, true-zero/underflow handling and partial-zero priors. Schema2 retains explicit no-CS and original component IDs. Locus plotting adds source-bound local annotation and coordinated LD heatmaps. See docs/upstream_differences.md.
