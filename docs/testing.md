# Tests and development

```sh
python -m pip install -c requirements-dev.lock '.[jit,annotation,test]'
python -m pytest tests
python examples/quickstart.py
python examples/check_example.py --output-dir example-output
python -m pip install -r tools/docs-requirements.txt
python tools/build_docs.py --check
python tools/check_docs.py
python -m pip install build
python -m build
```

Public tests include the frozen upstream synthetic 500-SNP fixture and hand-written analytical fixtures: independent high-precision
configuration enumeration, ABF equations, sdY regression, prior scaling/zeros,
SNP permutations/intersections, nonconsecutive CS components, null columns,
zero evidence and underflow, NumPy/Numba parity, runtime independence, plotting
semantics and local miniature GTF/GFF3 indexing. They need no private dataset,
network data download, R installation or original development workspace.

Every plotted probability is checked against its result object. Plotting tests
also exercise source/cache identity, coordinate boundaries, LD ordering and
SVG/PDF/PNG behavior. Optional annotation/JIT dependencies should be installed
for the full public suite. No external real dataset is silently skipped as a
passing scientific comparison.

Private release evidence separately consumes all100 frozen regional native and
official-R fits through the unchanged pycoloc inference implementation, including
no-CS and nonconvergence. It records component-pair additions/removals, H0–H4 and
SNP.PP.H4 changes and ranks. Fine-mapping acceptance is PIP atol1e-5/rtol0;
downstream diagnostic differences do not add an undeclared rejection threshold.
Earlier fixed-BF mathematical checks use tighter tolerances where appropriate.

The supplied GitHub workflow runs the public suite/build/docs/examples locally
reproducibly. It contains no publishing job or secrets. Local success does not
claim that GitHub Actions has executed.
