# Installation

Install this release candidate and its bundled example from source:

```sh
git clone https://github.com/LucaJiang/pycoloc.git
cd pycoloc
python -m venv .venv
. .venv/bin/activate
python -m pip install .
python examples/quickstart.py
```

The package requires Python ≥3.10, NumPy ≥1.24, SciPy ≥1.10, pandas ≥2 and
matplotlib ≥3.7. Executed checks cover CPython 3.12 on Linux x86-64. Keep the
checkout for examples and frozen reference data; installation can download
dependencies, then the examples run offline.

ABF, BF pairing and consuming existing SuSiE fits require only pycoloc's
Python dependencies. For native fine-mapping, install
[prusie](https://github.com/LucaJiang/prusie) separately and follow the
[fit-reuse example](example.md#fit-once-and-reuse-the-result).

## Optional features

From the activated environment and checkout:

```sh
python -m pip install '.[jit]'         # Numba hypothesis kernel
python -m pip install '.[annotation]'  # Local GTF/GFF3 indexing
```

Numba compilation contributes to the first call; measure it separately from
warm calls. Prepared annotation models render without gffutils. Validation and
package-building commands are in [contributing](contributing.md).
