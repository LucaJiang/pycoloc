# pycoloc 0.2.2rc3 (release candidate)

The current export uses **prusie** for the repository, distribution and import.
See the [rename and migration notes](migration.md). Historical performance and
R-agreement records below retain their original binary identities.

- Preserve existing corrected inference, locus plotting and local annotation behavior.
- Bundle the complete official synthetic 500-SNP teaching data, portable pinned-R fits, distinct mathematical/native expectations, and a checksum-validating offline public self-check.
- Add executed ABF/identical-fit pairing comparisons, explicit weighted counterexamples, no-CS and nonconvergence cases, and an optional native package-pair mode.
- Add a static `docs/` Pages site with Markdown sources, local assets, deterministic rendering and offline link/anchor checks.
- Extend local CI to run the frozen example, failure tests and rendered-doc freshness check.

The new version identifies documentation, examples and packaging changes. Statistical implementation and inference defaults are unchanged from pycoloc 0.2.2rc1. This package version does not imply that any prusie optimization is accepted. Its separately selected source must pass the majority-case runtime policy; B01 remains historical experimental evidence.

The source includes GitHub CI checks; remote CI success, a DOI, package-index publication and default-environment replacement are not claimed by the local validation. Local execution is tested on Linux x86-64 / CPython 3.12. Other platforms require supported dependency builds and are not claimed as locally tested.

The [agreement report](r_agreement.md) distinguishes the public synthetic example from the independently reviewed private real-region panel. Primary known-n fine-mapping reference settings retain max_iter=100 and PIP absolute error ≤1e-5, rtol=0, plus input/model/probability validity. Fixed-BF pycoloc checks retain the existing tighter thresholds. Corrected upstream differences and nonconvergence remain explicit limitations.

## GitHub source delivery

Repository: [pycoloc](https://github.com/LucaJiang/pycoloc). The publication changes add repository links and Pages setup instructions; numerical source, frozen inputs and expected outputs are unchanged.

The previously validated 0.2.2rc3 wheels and source distributions are **pre-publication snapshots**. They predate the repository metadata and documentation links in this checkout. Their original hashes and installation evidence are retained; this GitHub delivery does not create replacement distributions with the same filenames or claim a PyPI release. Install this checkout for the maintained metadata and documentation.
