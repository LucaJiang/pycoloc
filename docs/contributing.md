# Contributing

Start with the checkout and environment in [installation](install.md). Run the
existing tests for affected code, and document any changed scientific behavior.

## Tests and packages

```sh
python -m pip install '.[jit,annotation,test]' ruff==0.16.9
ruff check
ruff format --check
python -m pytest tests
python examples/quickstart.py
python examples/check_example.py --output-dir example-output
python -m pip install build
python -m build
```

The build writes a wheel and source distribution into `dist/`; use the actual
filename when installing your own wheel. The source distribution includes the
repository examples. These commands do not publish artifacts.

The public tests cover independent configuration enumeration, ABF equations,
weighted alignment, explicit priors, null components, partial zeros and stable
H3, NumPy/Numba parity, runtime independence, plotting and local miniature
GTF/GFF3 indexing. Install the JIT and annotation extras for the full suite.
`tests/test_example.py` includes checksum and numerical failure controls and
verifies standalone execution without R or prusie. Optional native integration
is documented in the [example guide](example.md#optional-native-fitting-with-prusie).

## Documentation

Markdown in `docs/` is the source for the committed static HTML. From the
repository root:

```sh
python -m pip install -r tools/docs-requirements.txt
python tools/build_docs.py
python tools/build_docs.py --check
python tools/check_docs.py
python -m http.server 8000 --bind 127.0.0.1 --directory docs
```

Open `http://127.0.0.1:8000/`, or open `docs/index.html` directly. Navigation,
code and tables work without JavaScript; Tab reaches the skip link and scrollable
content. Links are relative so the site also works below a repository subpath.
`--check` detects stale HTML; the link checker validates local pages, anchors,
and assets. Review desktop and mobile rendering after layout changes.
When merging a page, remove its obsolete HTML and update inbound links too.

The published site is [pycoloc documentation](https://lucajiang.github.io/pycoloc/).
An owner can configure it under **Settings → Pages → Build and deployment**:
**Deploy from a branch → main → /docs → Save**. `docs/.nojekyll` serves the
prebuilt files. Committing source and deploying Pages are separate operations;
check deployment status after publishing changes.

## Scientific checks and provenance

The [R agreement report](r_agreement.md) identifies the measured versions,
parameters, tolerances, corrected behavior and limits. Keep official-R outputs,
mathematical expectations and native regression snapshots distinct. Never
refresh frozen outputs merely to make a failed check pass. Reference regeneration
uses a separate directory and the pinned recipe in
[the data README](https://github.com/LucaJiang/pycoloc/blob/main/examples/data/README.md).
See [source mapping](reference_mapping.md) for implementation responsibilities.

Retain GPL and third-party notices with code and data. Cite the statistical
methods as well as the software version. Tests use redistributable analytical
inputs and the fixed synthetic toy example; they do not need private study
data or an R installation. The workflow in `.github/workflows/ci.yml` runs local
build/test/doc commands on GitHub; consult the actual workflow run for CI status.
