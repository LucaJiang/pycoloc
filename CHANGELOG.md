# Changelog

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
