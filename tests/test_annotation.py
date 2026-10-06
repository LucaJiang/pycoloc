"""Small local files exercise actual SQLite relationships, not mocked parsing."""

import json
import socket

import pytest

from pycoloc.annotation import (
    annotation_from_dict,
    annotation_to_dict,
    prepare_region_annotation,
    union_exons,
)
from pycoloc.annotation.models import select_transcript


@pytest.fixture
def annotation_source(tmp_path, monkeypatch):
    def offline(*args, **kwargs):
        raise AssertionError("Annotation attempted network access")

    monkeypatch.setattr(socket.socket, "connect", offline)
    path = tmp_path / "mini.gtf"
    rows = [
        ("gene", 50, 600, 'gene_id "g1"; gene_name "TARGET";'),
        ("transcript", 50, 600, 'gene_id "g1"; transcript_id "t1";'),
        ("exon", 50, 100, 'gene_id "g1"; transcript_id "t1";'),
        ("exon", 200, 250, 'gene_id "g1"; transcript_id "t1";'),
        ("transcript", 70, 550, 'gene_id "g1"; transcript_id "t2"; tag "MANE_Select";'),
        ("exon", 70, 110, 'gene_id "g1"; transcript_id "t2";'),
        ("exon", 251, 280, 'gene_id "g1"; transcript_id "t2";'),
        ("gene", 350, 500, 'gene_id "g2"; gene_name "CONTEXT";'),
        ("transcript", 350, 500, 'gene_id "g2"; transcript_id "t3";'),
        ("exon", 350, 400, 'gene_id "g2"; transcript_id "t3";'),
        ("gene", 800, 900, 'gene_id "g3"; gene_name "OUTSIDE";'),
    ]
    path.write_text(
        "##genome-build GRCh38\n"
        + "".join(f"chr22\ttest\t{f}\t{a}\t{b}\t.\t+\t.\t{attrs}\n" for f, a, b, attrs in rows)
    )
    options = dict(
        chromosome="22",
        start=90,
        end=450,
        genome_build="GRCh38",
        version="fixture-1",
        source_name="local test",
        cache_dir=tmp_path / "cache",
        index_scope=["chr22"],
    )
    return path, options


def test_gtf_boundary_relationships_and_union(annotation_source):
    path, options = annotation_source
    a = prepare_region_annotation(path, **options)
    assert [g.gene_id for g in a.genes] == ["g1", "g2"]
    g = a.genes[0]
    assert (g.start, g.end) == (50, 600)  # not clipped at the query boundaries
    assert [t.transcript_id for t in g.transcripts] == ["t1", "t2"]
    assert len(g.exons) == 4 and len(g.transcripts[0].exons) == 2
    assert union_exons(g.exons) == ((50, 110), (200, 280))
    transcript, rule = select_transcript(g)
    assert transcript.transcript_id == "t2" and rule == "annotated MANE_Select"
    assert annotation_from_dict(annotation_to_dict(a)) == a


def test_reuse_invalidation_and_corruption(annotation_source):
    path, options = annotation_source
    first = prepare_region_annotation(path, **options)
    from pathlib import Path

    db = Path(first.provenance["db_path"])
    before = db.stat().st_mtime_ns
    second = prepare_region_annotation(path, **options)
    assert second.provenance["cache_reused"] and db.stat().st_mtime_ns == before
    assert first.genes == second.genes
    meta = db.with_suffix(".json")
    value = json.loads(meta.read_text())
    value["identity"]["version"] = "stale"
    meta.write_text(json.dumps(value))
    repaired = prepare_region_annotation(path, **options)
    assert not repaired.provenance["cache_reused"] and repaired.genes == first.genes
    path.write_text(path.read_text() + "# Source changed\n")
    changed = prepare_region_annotation(path, **options)
    assert changed.provenance["cache_key"] != first.provenance["cache_key"]
    assert changed.provenance["source_sha256"] != first.provenance["source_sha256"]


def test_chromosome_scope_and_build(annotation_source):
    path, options = annotation_source
    for update, match in [
        ({"chromosome": "21"}, "scope"),
        ({"genome_build": "GRCh37"}, "GRCh38"),
        ({"chromosome_policy": "exact", "index_scope": None}, "Chromosome absent"),
    ]:
        with pytest.raises(ValueError, match=match):
            prepare_region_annotation(path, **{**options, **update})
    path.write_text(path.read_text().replace("GRCh38", "GRCh37"))
    with pytest.raises(ValueError, match="header genome build"):
        prepare_region_annotation(path, **options)


def test_gff3_parent_links(tmp_path, monkeypatch):
    monkeypatch.setattr(socket.socket, "connect", lambda *a, **k: pytest.fail("network attempted"))
    path = tmp_path / "mini.gff3"
    path.write_text(
        "##gff-version 3\n##genome-build GRCh38\n"
        "22\ttest\tgene\t50\t600\t.\t-\t.\tID=g;Name=Target\n"
        "22\ttest\tmRNA\t60\t590\t.\t-\t.\tID=t;Parent=g;tag=Ensembl_canonical\n"
        "22\ttest\texon\t60\t120\t.\t-\t.\tID=e1;Parent=t\n"
        "22\ttest\texon\t300\t400\t.\t-\t.\tID=e2;Parent=t\n"
    )
    a = prepare_region_annotation(
        path,
        chromosome="chr22",
        start=100,
        end=200,
        genome_build="GRCh38",
        version="fixture",
        cache_dir=tmp_path / "cache",
    )
    assert len(a.genes) == 1 and len(a.genes[0].transcripts) == 1
    assert [e.exon_id for e in a.genes[0].transcripts[0].exons] == ["e1", "e2"]
    assert len(a.genes[0].exons) == 2 and a.genes[0].strand == "-"
    assert select_transcript(a.genes[0])[1] == "annotated Ensembl_canonical"


def test_recognized_gencode_version_rejected(annotation_source):
    path, options = annotation_source
    path.write_text(
        "##description: human genome (GRCh38), version 50 (Ensembl 116)\n" + path.read_text()
    )
    with pytest.raises(ValueError, match="header version"):
        prepare_region_annotation(path, **{**options, "source_name": "GENCODE", "version": "49"})
    a = prepare_region_annotation(path, **{**options, "source_name": "GENCODE", "version": "50"})
    assert len(a.genes) == 2


def test_source_mutation_during_build_is_not_committed(annotation_source, monkeypatch):
    import gffutils

    path, options = annotation_source
    original = gffutils.create_db

    def mutate(*a, **k):
        result = original(*a, **k)
        path.write_text(path.read_text() + "# changed during indexing\n")
        return result

    monkeypatch.setattr(gffutils, "create_db", mutate)
    with pytest.raises(ValueError, match="changed during index"):
        prepare_region_annotation(path, **options)
    assert not list(options["cache_dir"].glob("*.sqlite"))
    assert not list(options["cache_dir"].glob("*.json"))


@pytest.mark.parametrize("rna_type", ["transcript", "ncRNA", "lnc_RNA"])
def test_gencode_gff3_rna_parent_types(tmp_path, rna_type):
    path = tmp_path / "gencode.gff3"
    path.write_text(
        "##gff-version 3\n##description: human genome (GRCh38), version 50 (Ensembl 116)\n"
        "chr22\tGENCODE\tgene\t50\t600\t.\t+\t.\tID=ENSG1;gene_name=TARGET;gene_type=lncRNA\n"
        f"chr22\tGENCODE\t{rna_type}\t60\t590\t.\t+\t.\tID=ENST1;Parent=ENSG1;tag=MANE_Select;transcript_name=TARGET-1\n"
        "chr22\tGENCODE\texon\t60\t120\t.\t+\t.\tID=e1;Parent=ENST1\n"
    )
    a = prepare_region_annotation(
        path,
        chromosome="22",
        start=100,
        end=200,
        genome_build="GRCh38",
        version="50",
        cache_dir=tmp_path / "cache",
    )
    assert a.genes[0].name == "TARGET" and a.genes[0].transcripts[0].name == "TARGET-1"
    assert a.genes[0].transcripts[0].exons[0].start == 60
    assert select_transcript(a.genes[0])[1] == "annotated MANE_Select"
