import builtins
import os
import subprocess

from pycoloc import coloc_abf, coloc_susie


def test_runtime_without_r_or_prusie(monkeypatch, datasets, fits):
    original = builtins.__import__

    def guarded(name, *args, **kwargs):
        if name.split(".")[0] in {"prusie", "rpy2"}:
            raise AssertionError("unexpected fine-mapping or R import")
        return original(name, *args, **kwargs)

    def denied(*args, **kwargs):
        raise AssertionError("statistical runtime attempted external process")

    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(builtins, "__import__", guarded)
    monkeypatch.setattr(subprocess, "Popen", denied)
    monkeypatch.setattr(os, "system", denied)
    assert coloc_abf(*datasets).nsnps == 4
    assert coloc_susie(*fits).nsnps == 4
