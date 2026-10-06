"""Offline gffutils SQLite indexing with content-addressed, locked provenance."""

import gzip
import hashlib
import importlib.metadata
import json
import os
import re
import sqlite3
import tempfile
from contextlib import closing
from dataclasses import asdict
from pathlib import Path

from .models import (
    ExonModel,
    GeneModel,
    RegionAnnotation,
    TranscriptModel,
    normalize_chromosome,
    validate_region,
)

INDEX_SCHEMA = "pycoloc-annotation-2"
TRANSCRIPT_TYPES = {
    "transcript",
    "mRNA",
    "ncRNA",
    "lnc_RNA",
    "tRNA",
    "rRNA",
    "miRNA",
    "snRNA",
    "snoRNA",
    "scRNA",
    "primary_transcript",
    "pseudogenic_transcript",
    "processed_transcript",
}


def _hash(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _header_metadata(path):
    opener = gzip.open if path.suffix == ".gz" else open
    builds, versions = set(), set()
    with opener(path, "rt") as stream:
        for _ in range(200):
            line = stream.readline()
            if not line.startswith("#"):
                break
            builds.update(re.findall(r"GRCh\d+", line))
            if "description" in line:
                versions.update(re.findall(r"\bversion\s+(\d+)\b", line))
    return builds, versions


def prepare_region_annotation(
    source,
    *,
    chromosome,
    start,
    end,
    genome_build,
    version,
    cache_dir,
    source_name="GENCODE",
    chromosome_policy="chr",
    index_scope=None,
    source_provenance=None,
):
    """Index a local GTF/GFF3 and query complete genes overlapping a region.

    No network access is performed. Install ``pycoloc[annotation]`` to index;
    prepared RegionAnnotation objects need no gffutils to plot. ``index_scope``
    lists all permitted chromosomes for a deliberate chromosome subset, otherwise
    the file is declared a complete supplied source. Coordinate-window subsets
    are not supported. ``source_provenance`` records an upstream source/hash and
    deterministic subsetting rule without replacing the actual indexed-file hash.
    A changed source/version/build/indexer/schema gets a separate index identity.
    Corrupt or mismatched cache files are rebuilt while holding a POSIX file lock.
    """
    validate_region(chromosome, start, end, genome_build)
    chromosome = normalize_chromosome(chromosome, chromosome_policy)
    scope = (
        None
        if index_scope is None
        else sorted({normalize_chromosome(c, chromosome_policy) for c in index_scope})
    )
    if scope is not None and chromosome not in scope:
        raise ValueError("Requested chromosome is outside the declared index scope")
    source = Path(source).expanduser().resolve(strict=True)
    builds, header_versions = _header_metadata(source)
    if builds and builds != {genome_build}:
        raise ValueError("Annotation header genome build disagrees with requested build")
    if source_name.upper() == "GENCODE" and header_versions and header_versions != {str(version)}:
        raise ValueError("GENCODE header version disagrees with supplied version")
    if not str(version) or not str(source_name):
        raise ValueError("Explicit annotation source and version are required")
    name = source.name.removesuffix(".gz").lower()
    if name.endswith(".gtf"):
        dialect = "gtf"
    elif name.endswith((".gff", ".gff3")):
        dialect = "gff3"
    else:
        raise ValueError("Local annotation must be .gtf, .gff or .gff3, optionally gzipped")
    try:
        import gffutils
    except ImportError as exc:
        raise ImportError('Indexing requires pip install "pycoloc[annotation]"') from exc
    import fcntl

    identity = dict(
        source=source_name,
        version=str(version),
        genome_build=genome_build,
        source_path=str(source),
        source_sha256=_hash(source),
        format=dialect,
        indexer="gffutils",
        indexer_version=importlib.metadata.version("gffutils"),
        schema=INDEX_SCHEMA,
        chromosome_policy=chromosome_policy,
        index_scope=scope,
        upstream=source_provenance or {},
        coordinate_convention="1-based inclusive",
        inference="disabled; complete explicit gene/transcript records required",
    )
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    cache = Path(cache_dir).expanduser().resolve()
    cache.mkdir(parents=True, exist_ok=True)
    dbpath, metapath = cache / (key + ".sqlite"), cache / (key + ".json")
    reused = False
    with (cache / (key + ".lock")).open("a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if dbpath.is_file() and metapath.is_file():
            try:
                previous = json.loads(metapath.read_text())
                with closing(sqlite3.connect(f"file:{dbpath}?mode=ro", uri=True)) as connection:
                    healthy = connection.execute("PRAGMA quick_check").fetchone()[0] == "ok"
                reused = (
                    healthy
                    and previous["identity"] == identity
                    and previous["db_sha256"] == _hash(dbpath)
                )
            except (ValueError, KeyError, sqlite3.DatabaseError):
                reused = False
        if not reused:
            fd, tempname = tempfile.mkstemp(prefix=key + ".", suffix=".sqlite", dir=cache)
            os.close(fd)
            tmp = Path(tempname)
            tmpmeta = tmp.with_suffix(".json")
            db = None
            try:
                db = gffutils.create_db(
                    str(source),
                    dbfn=str(tmp),
                    force=True,
                    disable_infer_genes=True,
                    disable_infer_transcripts=True,
                    merge_strategy="error",
                    keep_order=True,
                    pragmas={
                        "journal_mode": "DELETE",
                        "synchronous": "FULL",
                        "cache_size": 10000,
                        "page_size": 4096,
                    },
                )
                seqids = sorted(set(db.seqids()))
                if scope is not None and {
                    normalize_chromosome(c, chromosome_policy) for c in seqids
                } != set(scope):
                    raise ValueError("Declared index scope differs from actual source chromosomes")
                if db.count_features_of_type("gene") == 0:
                    raise ValueError("Annotation source has no explicit gene records")
                db.conn.close()
                if _hash(source) != identity["source_sha256"]:
                    raise ValueError(
                        "Annotation source changed during index creation; retry with a stable source"
                    )
                previous = dict(
                    identity=identity, db_sha256=_hash(tmp), seqids=seqids, cache_key=key
                )
                tmpmeta.write_text(json.dumps(previous, indent=2, sort_keys=True) + "\n")
                os.replace(tmp, dbpath)
                os.replace(tmpmeta, metapath)
            finally:
                if db is not None:
                    db.conn.close()
                tmp.unlink(missing_ok=True)
                tmpmeta.unlink(missing_ok=True)
        # Read only: FeatureDB ordinarily accepts a read-only sqlite connection.
        connection = sqlite3.connect(f"file:{dbpath}?mode=ro", uri=True)
        try:
            db = gffutils.FeatureDB(connection, keep_order=True)
            candidates = [
                c
                for c in previous["seqids"]
                if normalize_chromosome(c, chromosome_policy) == chromosome
            ]
            if len(candidates) != 1:
                raise ValueError("Chromosome absent from index or normalization is ambiguous")
            seqid = candidates[0]

            def attr(feature, *keys):
                for k in keys:
                    values = feature.attributes.get(k, [])
                    if values:
                        return values[0]
                return None

            def exon(feature):
                return ExonModel(
                    feature.id, chromosome, int(feature.start), int(feature.end), feature.strand
                )

            genes = []
            for gene in db.region(
                seqid=seqid,
                start=int(start),
                end=int(end),
                featuretype="gene",
                completely_within=False,
            ):
                transcripts = []
                for t in db.children(gene, level=1, order_by=("start", "end")):
                    if t.featuretype not in TRANSCRIPT_TYPES:
                        if any(db.children(t, featuretype="exon")):
                            raise ValueError(
                                "Unsupported exon parent feature type: " + t.featuretype
                            )
                        continue
                    exons = tuple(
                        exon(e)
                        for e in db.children(t, featuretype="exon", order_by=("start", "end"))
                    )
                    transcripts.append(
                        TranscriptModel(
                            t.id,
                            chromosome,
                            t.start,
                            t.end,
                            t.strand,
                            exons,
                            tuple(t.attributes.get("tag", [])),
                            attr(t, "transcript_name", "Name"),
                            attr(t, "transcript_type", "transcript_biotype", "biotype"),
                        )
                    )
                exons = tuple(
                    exon(e)
                    for e in db.children(gene, featuretype="exon", order_by=("start", "end"))
                )
                genes.append(
                    GeneModel(
                        gene.id,
                        attr(gene, "gene_name", "Name") or gene.id,
                        chromosome,
                        gene.start,
                        gene.end,
                        gene.strand,
                        tuple(transcripts),
                        exons,
                        attr(gene, "gene_type", "gene_biotype", "biotype"),
                    )
                )
        finally:
            connection.close()
    provenance = dict(
        identity,
        db_path=str(dbpath),
        db_sha256=previous["db_sha256"],
        cache_key=key,
        cache_reused=reused,
        actual_seqids=previous["seqids"],
    )
    return RegionAnnotation(
        chromosome,
        int(start),
        int(end),
        genome_build,
        tuple(sorted(genes, key=lambda g: (g.start, g.end, g.gene_id))),
        provenance,
        chromosome_policy,
    )


def annotation_to_dict(annotation):
    """Return a JSON-serializable complete annotation snapshot."""
    return asdict(annotation)


def annotation_from_dict(value):
    """Load a prepared snapshot without gffutils, files, or network requests."""
    genes = []
    for g in value["genes"]:
        transcripts = tuple(
            TranscriptModel(
                **{
                    **t,
                    "exons": tuple(ExonModel(**e) for e in t["exons"]),
                    "tags": tuple(t["tags"]),
                }
            )
            for t in g["transcripts"]
        )
        genes.append(
            GeneModel(
                **{
                    **g,
                    "transcripts": transcripts,
                    "exons": tuple(ExonModel(**e) for e in g["exons"]),
                }
            )
        )
    return RegionAnnotation(**{**value, "genes": tuple(genes)})
