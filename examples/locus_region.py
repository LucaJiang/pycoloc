"""Render a local locus bundle with installed pycoloc, without fitting or downloads.

Usage: python locus_region.py --config region1_susie_bundle.json --output-prefix figures/locus
The JSON bundle schema and direct Python API are documented in docs/locus_plot.md.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def render_bundle(config_path, output_prefix, *, max_heatmap_snps=None, fontsize=7.5):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycoloc import LDMatrix, coloc_abf, coloc_susie, plot_coloc_locus
    from pycoloc.annotation import annotation_from_dict, prepare_region_annotation

    path = Path(config_path).resolve()
    config = json.loads(path.read_text())

    def resolve(value):
        p = Path(value)
        return p if p.is_absolute() else path.parent / p

    datasets = [json.loads(resolve(p).read_text()) for p in config["datasets"]]
    options = config.get("coloc_options", {})
    if config["method"] == "abf":
        result = coloc_abf(*datasets, **options)
    elif config["method"] == "susie":
        fits = []
        for entry in config["fits"]:
            with np.load(resolve(entry["npz"]), allow_pickle=False) as arrays:
                fits.append({**json.loads(resolve(entry["json"]).read_text()), **dict(arrays)})
        result = coloc_susie(*fits, **options)
    else:
        raise ValueError("Bundle method must be abf or susie")
    region = config["region"]
    if "annotation_snapshot" in config:
        annotation = annotation_from_dict(
            json.loads(resolve(config["annotation_snapshot"]).read_text())
        )
    else:
        source = config["annotation_source"]
        annotation = prepare_region_annotation(
            resolve(source["path"]),
            **region,
            version=source["version"],
            cache_dir=resolve(source["cache_dir"]),
            source_name=source.get("source_name", "GENCODE"),
            index_scope=source.get("index_scope"),
            source_provenance=source.get("provenance"),
        )
    panels = {}
    for name in ("gwas", "eqtl"):
        spec = config.get("ld", {}).get(name)
        if spec is None:
            panels[name] = None
            continue
        with np.load(resolve(spec["npz"]), allow_pickle=False) as arrays:
            panels[name] = LDMatrix(
                arrays[spec["values_key"]],
                arrays[spec["snp_ids_key"]],
                arrays[spec["positions_key"]],
                kind=spec["kind"],
                panel=spec["panel"],
                provenance=spec.get("provenance", {}),
            )
    heatmap = config.get("ld", {}).get("heatmap")
    if heatmap not in (None, "gwas", "eqtl"):
        raise ValueError("Bundle heatmap must explicitly select gwas, eqtl or null")
    plot_options = config.get("plot_options", {}).copy()
    if max_heatmap_snps is not None:
        plot_options["max_heatmap_snps"] = max_heatmap_snps
    fig, axes = plot_coloc_locus(
        result,
        *datasets,
        region=(region["chromosome"], region["start"], region["end"], region["genome_build"]),
        annotation=annotation,
        target_gene=config.get("target_gene"),
        signal_pair=config.get("signal_pair"),
        leads=config.get("leads"),
        lead_rule=config.get("lead_rule", "min_p"),
        trait1=config.get("trait1", "GWAS"),
        trait2=config.get("trait2", "eQTL"),
        gwas_ld=panels["gwas"],
        eqtl_ld=panels["eqtl"],
        heatmap_ld=panels.get(heatmap),
        output_prefix=output_prefix,
        fontsize=fontsize,
        **plot_options,
    )
    metadata = fig.coloc_metadata
    plt.close(fig)
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Local locus bundle JSON")
    parser.add_argument("--output-prefix", required=True, help="SVG/PDF/PNG/JSON output stem")
    parser.add_argument("--max-heatmap-snps", type=int, default=None)
    parser.add_argument("--fontsize", type=float, default=7.5)
    args = parser.parse_args()
    metadata = render_bundle(
        args.config,
        args.output_prefix,
        max_heatmap_snps=args.max_heatmap_snps,
        fontsize=args.fontsize,
    )
    print(json.dumps({"outputs": metadata["outputs"], "extent_pass": metadata["extent_pass"]}))


if __name__ == "__main__":
    main()
