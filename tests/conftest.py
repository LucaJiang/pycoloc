import numpy as np
import pytest


@pytest.fixture
def datasets():
    common = dict(
        snp=["a", "b", "c", "d"],
        type="cc",
        position=[100, 200, 300, 400],
        chromosome="22",
        N=1000,
        s=0.25,
        MAF=[0.1, 0.2, 0.3, 0.4],
    )
    d1 = dict(
        common,
        beta=[0.1, 0.4, -0.3, 0.2],
        varbeta=[0.02, 0.015, 0.04, 0.025],
        pvalues=[0.47, 0.001, 0.13, 0.2],
    )
    d2 = dict(
        common,
        beta=[0.2, 0.5, 0.1, -0.1],
        varbeta=[0.03, 0.02, 0.02, 0.04],
        pvalues=[0.25, 0.0004, 0.48, 0.62],
    )
    return d1, d2


@pytest.fixture
def fits():
    # Deliberate nonconsecutive CS indices; row 1 is not a detected signal.
    f1 = dict(
        variant_ids=["a", "b", "c", "d"],
        lbf_variable=np.array([[1.0, 4.0, 0.0, 2.0], [9.0, 0.0, 0.0, 0.0], [0.0, 1.0, 5.0, 2.0]]),
        converged=True,
        sets=dict(cs=[[1], [2]], cs_index=[0, 2]),
        null_index=None,
    )
    f2 = dict(
        variant_ids=["a", "b", "c", "d"],
        lbf_variable=np.array([[0.0, 3.0, 1.0, 2.0], [0.0, 2.0, 6.0, 0.0]]),
        converged=True,
        sets=dict(cs=[[1], [2]], cs_index=[0, 1]),
        null_index=None,
    )
    return f1, f2
