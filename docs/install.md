# Installation

Install pycoloc **0.2.2** with Python ≥3.10 from the GitHub Release wheel:

```sh
python -m pip install https://github.com/LucaJiang/pycoloc/releases/download/v0.2.2/pycoloc-0.2.2-py3-none-any.whl
```

Or install the tagged source:

```sh
python -m pip install "git+https://github.com/LucaJiang/pycoloc.git@v0.2.2"
```

pycoloc is not published on PyPI. The PyPI project named `pycoloc` is unrelated:
`pip install pycoloc` installs that other package. This GitHub repository is private;
access to its source and Release assets requires authorization. For the wheel path,
download the asset while signed in and install the local file if pip cannot authenticate.

## Build from source and run the examples

```sh
git clone --branch v0.2.2 https://github.com/LucaJiang/pycoloc.git
cd pycoloc
python -m venv .venv
. .venv/bin/activate
python -m pip install .
python examples/quickstart.py
```

The package requires Python ≥3.10, NumPy ≥1.24, SciPy ≥1.10, pandas ≥2 and
matplotlib ≥3.7. Keep the checkout for examples and frozen reference data.
Once dependencies are installed, the examples run offline.

## Optional features

For native fine-mapping, install prusie from PyPI (the integration is checked with 0.3.1):

```sh
python -m pip install prusie
```

ABF and pairing existing fits do not require prusie. See the
[fit-reuse example](example.md#fit-once-and-reuse-the-result).
From a checkout, install optional Numba or annotation support with:

```sh
python -m pip install '.[jit]'
python -m pip install '.[annotation]'
```

Numba compiles on first use. Prepared annotation models render without gffutils.
See [contributing](contributing.md) for tests and package builds.
