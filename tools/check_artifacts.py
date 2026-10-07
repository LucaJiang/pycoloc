#!/usr/bin/env python3
"""Inspect actual wheel/sdist metadata and license payloads without installing."""
from __future__ import annotations

import argparse
from email.parser import Parser
import hashlib
import json
from pathlib import Path
import tarfile
import zipfile


def inspect(path):
    path = Path(path)
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            metadata = Parser().parsestr(archive.read(next(n for n in names if n.endswith(".dist-info/METADATA"))).decode())
            wheel = Parser().parsestr(archive.read(next(n for n in names if n.endswith(".dist-info/WHEEL"))).decode())
            assert wheel.get("Root-Is-Purelib") == "true", "Wheel is not pure Python"
            assert wheel.get_all("Tag") == ["py3-none-any"], "Unexpected wheel compatibility tags"
            assert path.name.endswith("-py3-none-any.whl"), "Wheel filename disagrees with tag"
            assert any(n == "pycoloc/abf.py" for n in names), "Missing installed public module"
            required = ["LICENSE", "THIRD_PARTY_NOTICES", "examples/data/LICENSE.GPL-3",
                        "examples/data/upstream/coloc_DESCRIPTION"]
            prefix = next(n.split("/licenses/")[0] + "/licenses/" for n in names if "/licenses/" in n)
            for item in required:
                assert prefix + item in names, f"Wheel missing license payload: {item}"
            assert b"GNU GENERAL PUBLIC LICENSE" in archive.read(prefix + "LICENSE")
            assert b"GNU GENERAL PUBLIC LICENSE" in archive.read(prefix + "examples/data/LICENSE.GPL-3")
            assert b"Chris Wallace" in archive.read(prefix + "THIRD_PARTY_NOTICES")
            assert b"6.0.3" in archive.read(prefix + "examples/data/upstream/coloc_DESCRIPTION")
            tags = wheel.get_all("Tag")
    else:
        with tarfile.open(path, "r:gz") as archive:
            names = archive.getnames()
            root = names[0].split("/")[0]
            required = ["pyproject.toml", "MANIFEST.in", "LICENSE", "THIRD_PARTY_NOTICES",
                        "examples/data/LICENSE.GPL-3", "examples/data/upstream/coloc_DESCRIPTION",
                        "examples/data/inputs.npz", "examples/check_example.py", "tests/test_correctness.py",
                        "benchmarks/coloc_case.py", "benchmarks/coloc_case.R", "benchmarks/monotonic_clock.c",
                        "tools/check_artifacts.py"]
            for item in required:
                assert root + "/" + item in names, f"sdist missing source/reference/license: {item}"
            metadata = Parser().parsestr(archive.extractfile(root + "/PKG-INFO").read().decode())
            tags = None
    assert metadata["Name"] == "pycoloc"
    assert metadata["License-Expression"] == "GPL-3.0-or-later"
    return {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size, "name": metadata["Name"], "version": metadata["Version"],
            "tags": tags, "license_expression": metadata["License-Expression"],
            "license_files": metadata.get_all("License-File"), "requires_dist": metadata.get_all("Requires-Dist"),
            "required_payloads_verified": required, "passed": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = [inspect(path) for path in args.artifacts]
    text = json.dumps(records, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text)
