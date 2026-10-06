# Installation

Install from the [GitHub source](https://github.com/LucaJiang/pycoloc). No PyPI
release or public prebuilt wheel is provided. The package declares Python ≥3.10
with NumPy ≥1.24, SciPy ≥1.10, pandas ≥2 and matplotlib ≥3.7. Executed checks
cover CPython 3.12 on Linux x86-64; other platforms need compatible dependencies
and their own validation. pycoloc itself is pure Python.

```sh
git clone https://github.com/LucaJiang/pycoloc.git
cd pycoloc
python -m venv .venv
. .venv/bin/activate
python -m pip install .
python examples/check_example.py --output-dir example-output
```

Keep the checkout to run the [bundled example](example.md). Installation can
download dependencies; the check runs offline afterward. R and prusie are not
needed for ABF, BF pairing, consuming existing fits, or plotting.

## Optional features

From the activated environment and checkout:

```sh
python -m pip install '.[jit]'         # Numba hypothesis kernel
python -m pip install '.[annotation]'  # Local GTF/GFF3 indexing with gffutils
```

Numba's first compilation takes additional time. Prepared annotation models
render without the annotation extra. For native fine-mapping, install
[prusie](https://github.com/LucaJiang/prusie) separately and follow the
[optional example](example.md#optional-native-fitting-with-prusie).
Tests and package-building commands are in [contributing](contributing.md).
