"""Analyze the complete bundled 500-SNP synthetic teaching example."""
import json
from pathlib import Path

import numpy as np

from pycoloc import coloc_abf

data = Path(__file__).resolve().parent / "data"
metadata = json.loads((data / "metadata.json").read_text())
with np.load(data / "inputs.npz", allow_pickle=False) as arrays:
    datasets = [
        dict(snp=metadata[key]["snp"], type=metadata[key]["type"],
             N=metadata[key]["N"], sdY=metadata[key]["sdY"],
             beta=arrays[key + "_beta"], varbeta=arrays[key + "_varbeta"])
        for key in ("D1", "D2")
    ]

result = coloc_abf(*datasets)
print(result.summary[[f"PP.H{i}.abf" for i in range(5)]])
print(result.results[["snp", "SNP.PP.H4"]].nlargest(5, "SNP.PP.H4"))
