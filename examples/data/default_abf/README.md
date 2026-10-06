# Literal default-prior ABF reference addendum

These two executed official `coloc::coloc.abf(D1, D2)` and
`coloc::coloc.abf(D3, D4)` calls supply no prior arguments. Actual coloc 6.0.3
defaults are p1=p2=0.0001 and p12=0.00001. This addendum supplements the
existing explicit-p12=0.000005 teaching comparisons without changing them.
All 500 synthetic upstream variables are retained. No SuSiE fit, simulation,
private-data extraction or timing was performed for this addendum.

The exact R, coloc, susieR, dependency, fixture and archive identities are in
the accompanying lock/environment files. The fixture is the same GPL-licensed
coloc teaching data documented in the shared example data README; that bundle
also supplies the original source/license and fixture checksum verifier.
R outputs retain jsonlite's documented decimal precision.

After acquiring and verifying the pinned sources and fixture, configure R_LIBS
to the frozen reference environment and run:

```sh
Rscript --vanilla regenerate_default_abf.R coloc_test_data.rda regenerated-defaults
```

Regeneration requires the pinned R environment; running the Python example with
the frozen JSON requires neither R nor network access. `checksums.json` binds
every public addendum file. Original machine-specific session information stays
in private task evidence and is not part of the distributable addendum.
