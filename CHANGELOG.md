# Changelog

## 0.2.2 — 2026-11-04

- First final release. The inference code is unchanged from 0.2.2rc3.
- Tested with prusie 0.3.1: its fit results go directly into `coloc_susie`,
  including refined, nonconverged and no-credible-set fits.
- New comparison with R: pycoloc after prusie 0.3.1 fits versus coloc after
  susieR 0.16.6 fits, on 100 real regions and the toy example.
- Distributed as a wheel and sdist on GitHub Releases. pycoloc is not on PyPI;
  the PyPI package named `pycoloc` is a different project.
- CI tests the installed wheel on Python 3.10 and 3.13, runs Ruff and builds the
  documentation. The third-party notice is now `THIRD_PARTY_NOTICES.md`.

## 0.2.2rc3

- Use prusie for the optional native-fitting example; standalone inference is unchanged.
- Clarify installation, examples and scientific validation in the documentation.

## 0.2.2rc2

- Add the offline 500-SNP synthetic teaching example with separate R, mathematical
  and native expectations, and an executed [agreement report](docs/r_agreement.md).
- Add static documentation and example integrity/failure checks.

## 0.2.2rc1

- Package the colocalisation, local annotation and locus plotting interfaces with
  self-contained tests and CI configuration.

## Earlier changes

- Correct weighted `SNP.PP.H4`, SNP-bound prior alignment and explicit effect-prior
  SD propagation. Support partial-zero priors and stable H3/zero/underflow handling.
- Result schema 2 preserves original component IDs and explicit unavailable states.
- Add locus plots with local annotation and coordinated LD heatmaps.

See [upstream differences](docs/upstream_differences.md) for the scientific
meaning and tests of the corrected fields.
