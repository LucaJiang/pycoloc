# Installation

pycoloc 0.2.2 requires Python 3.10 or newer, NumPy ≥1.24, SciPy ≥1.10, pandas ≥2
and matplotlib ≥3.7. It is pure Python, so one wheel works on every platform.

## From the GitHub Release

```sh
pip install https://github.com/LucaJiang/pycoloc/releases/download/v0.2.2/pycoloc-0.2.2-py3-none-any.whl
```

or from the tagged source:

```sh
pip install "git+https://github.com/LucaJiang/pycoloc.git@v0.2.2"
```

pycoloc is not on PyPI. The PyPI package named `pycoloc` is a different project,
so `pip install pycoloc` will not install this one.

## From a clone

Clone the repository if you want the examples and the stored R results:

```sh
git clone --branch v0.2.2 https://github.com/LucaJiang/pycoloc.git
cd pycoloc
python -m venv .venv
. .venv/bin/activate
pip install .
python examples/quickstart.py
```

The examples run offline.

## Optional dependencies

| Purpose | Command |
| --- | --- |
| Fit SuSiE models in Python (prusie 0.3.1 is tested) | `pip install prusie` |
| Numba kernel for the H3 sum | `pip install '.[jit]'` |
| Build gene annotation from GFF/GTF files with gffutils | `pip install '.[annotation]'` |

ABF and colocalisation of existing fits do not need prusie. Numba compiles on the
first call. Annotation models that are already prepared render without gffutils.
See [contributing](contributing.md) for running the tests and building packages.
