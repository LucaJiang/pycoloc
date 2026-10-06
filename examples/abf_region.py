"""Small hand-fixture API example; not a simulated or biological study."""

import matplotlib

matplotlib.use("Agg")
from pycoloc import coloc_abf, plot_coloc

d1 = dict(
    snp=["a", "b", "c"],
    position=[100, 200, 300],
    chromosome="22",
    beta=[0.1, 0.5, -0.2],
    varbeta=[0.01, 0.02, 0.03],
    type="cc",
)
d2 = dict(d1, beta=[0.0, 0.4, 0.1], type="quant", sdY=1.2)
result = coloc_abf(d1, d2, prior_weights1={"a": 1, "b": 2, "c": 1})
print(result.summary.to_string(index=False))
fig, axes = plot_coloc(result, d1, d2)
fig.savefig("hand_fixture.svg")
