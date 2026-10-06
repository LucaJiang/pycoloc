# Integration with prusie

The separate fine-mapping package now uses **prusie** for its repository,
distribution and Python import. **pycoloc 0.2.2rc3** updates the optional example,
import guards and documentation to that name. Its statistical implementation,
result schema, defaults and numerical tolerances are unchanged from 0.2.2rc2.

After installing the supplied local packages or their source trees:

```python
import pycoloc
import prusie
from importlib.metadata import version

assert pycoloc.__version__ == version("pycoloc")
assert prusie.__version__ == version("prusie")
```

Run from the pycoloc repository or extracted source distribution:

```sh
python examples/check_example.py --output-dir example-results
python examples/check_example.py --native-pair --output-dir pair-results
```

The default check and ordinary pycoloc computations do not require prusie or R.
The explicit `--native-pair` option imports `prusie` and reports
`prusie_version`. Replace old `import pyrsusie` statements in your own code;
the renamed package supplies no legacy import aliases. Its native extension is
`prusie._native`. Use `PRUSIE_NUM_THREADS` and `PRUSIE_SER_MATH` for native
runtime settings. No package-index publication is implied by these instructions.

The dataclass path changed from `pyrsusie.result.SusieResult` to
`prusie.result.SusieResult`; old pickles require their original environment.
Portable arrays and metadata remain compatible, and pycoloc still accepts
existing fit mappings without importing either fine-mapping package.

Frozen R and native expected outputs, fixture checksums and historical
performance facts retain their original bytes. The R-agreement report and
downloadable prior check records continue to identify pycoloc 0.2.2rc2 and
pyrsusie 0.2.3rc4 as the binaries actually used. Rename verification covers the
new installs, full test suites, offline examples, native pair, docs and source
equivalence; it is not a fresh independent statistical sign-off or benchmark.

## License status

The user accepted the existing **GPL-3.0-or-later** package licensing after
initially requesting MIT. Existing notices and upstream license/version terms
remain intact; this rename does not relicense third-party code or data.
The concrete retained provenance includes the native R 4.4 `Brent_fmin`
adaptation (GPL-2-or-later), pycoloc's adaptation from GPL coloc R files, and
the GPL coloc teaching fixture. Mathematical-method agreement alone is not
the reason for those notices. The earlier MIT preference is resolved by the
explicit GPL decision; no licensing rewrite is part of this release.
