# Installation

Python ≥3.10 is required. Runtime dependencies are NumPy≥1.24, SciPy≥1.10,
pandas≥2 and matplotlib≥3.7. Local release checks use CPython3.12/Linux x86-64.
The package is pure Python; dependency wheels must support your platform. Other
platforms are not claimed as locally tested.

From this repository root:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install .
python -m pip install -c requirements-dev.lock '.[jit,annotation,test]'
python -m pytest tests
python examples/quickstart.py
```

`jit` adds Numba for the optional hypothesis kernel; first compilation and warm
execution differ. `annotation` adds gffutils for local GTF/GFF3 indexing.
Prepared annotation models render without that extra. `test` adds pytest.
R and prusie are not runtime dependencies for ABF, BF pairing, consuming
existing fits, or plotting. The separate prusie repository builds distribution
prusie when native fine-mapping is wanted.

```sh
python -m pip install build
python -m build
python -m pip install --no-deps dist/pycoloc-*.whl
```

To install a supplied source archive, use its real filename with
`python -m pip install path/to/pycoloc-VERSION.tar.gz`. Building uses setuptools≥77.
No package-index release or remote CI success is implied by these local commands.

To obtain the source and its frozen offline example:

```sh
git clone https://github.com/LucaJiang/pycoloc.git
cd pycoloc
python -m pip install .
```

This installs from the checkout; it does not require or imply a PyPI release.
