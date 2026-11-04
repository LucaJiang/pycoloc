"""Public self-check must detect broken input/reference bundles and avoid hidden runtimes."""
import builtins
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / 'examples/check_example.py'
DATA = ROOT / 'examples/data'


def run_check(data, output):
    env = dict(os.environ, PATH='', PYTHONDONTWRITEBYTECODE='1')
    return subprocess.run([sys.executable, str(RUNNER), '--data-dir', str(data),
                           '--output-dir', str(output)], capture_output=True, text=True, env=env)


def test_self_check_cli_offline(tmp_path):
    result = run_check(DATA, tmp_path / 'result')
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((tmp_path / 'result/report.json').read_text())
    assert report['passed'] and report['integrity_files_checked'] >= 10
    assert report['fixture_snp_count'] == 500
    assert report['native_pair'] is False and report['prusie_version'] is None
    assert report['diagnostics']  # Raw upstream differences stay visible.


def test_tampered_input_checksum_fails(tmp_path):
    data = tmp_path / 'data'
    shutil.copytree(DATA, data)
    path = data / 'metadata.json'
    metadata = json.loads(path.read_text())
    metadata['D1']['sdY'] *= 2
    path.write_text(json.dumps(metadata))
    result = run_check(data, tmp_path / 'result')
    assert result.returncode != 0
    assert 'Checksum mismatch for metadata.json' in result.stdout
    assert not json.loads((tmp_path / 'result/report.json').read_text())['passed']


def test_tampered_expected_value_fails_numerically(tmp_path):
    data = tmp_path / 'data'
    shutil.copytree(DATA, data)
    path = data / 'r_coloc.json'
    reference = json.loads(path.read_text())
    reference['cases']['abf_D1_D2']['result']['summary']['PP.H4.abf'] -= 0.1
    path.write_text(json.dumps(reference))
    # Deliberately repair this one checksum to exercise the numerical check too.
    manifest = json.loads((data / 'checksums.json').read_text())
    manifest['files']['r_coloc.json']['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    (data / 'checksums.json').write_text(json.dumps(manifest))
    result = run_check(data, tmp_path / 'result')
    assert result.returncode != 0
    assert 'FAIL abf_D1_D2 / H0–H4' in result.stdout
    record = json.loads((tmp_path / 'result/report.json').read_text())
    assert not record['passed'] and 'error' not in record


def test_standalone_never_imports_finemapping(monkeypatch):
    spec = importlib.util.spec_from_file_location('public_example', RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original = builtins.__import__

    def guarded(name, *args, **kwargs):
        if name.split('.')[0] in {'prusie', 'rpy2'}:
            raise AssertionError('standalone example attempted an optional/R import')
        return original(name, *args, **kwargs)

    def denied(*args, **kwargs):
        raise AssertionError('standalone numerical checks attempted an external process')

    monkeypatch.setattr(builtins, '__import__', guarded)
    monkeypatch.setattr(subprocess, 'Popen', denied)
    monkeypatch.setenv('PATH', '')
    record = module.evaluate(DATA, expected_native=module.read_json(DATA / 'expected_native.json'))
    assert record['passed']


def test_version_string_alone_does_not_fail(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location('public_example_version_check', RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module.pycoloc, '__version__', 'future-version-label')
    assert module.main(['--output-dir', str(tmp_path / 'result')]) == 0
    report = json.loads((tmp_path / 'result/report.json').read_text())
    assert report['passed'] and report['version_differs']
