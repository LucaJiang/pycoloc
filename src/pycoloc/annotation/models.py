"""Provider-independent annotation models. All coordinates are 1-based inclusive."""

from dataclasses import dataclass, field
from math import isfinite
from numbers import Real
from typing import Any


def normalize_chromosome(value: str, policy: str = "chr") -> str:
    """Explicitly add/remove a leading 'chr', or require exact spelling."""
    value = str(value)
    if not value or policy not in {"chr", "none", "exact"}:
        raise ValueError("chromosome must be nonempty; policy is chr, none or exact")
    if policy == "exact":
        return value
    core = value[3:] if value.startswith("chr") else value
    return ("chr" if policy == "chr" else "") + core


def validate_region(chromosome, start, end, genome_build):
    if genome_build != "GRCh38":
        raise ValueError("Only explicitly identified Human GRCh38 is supported")
    _coordinates(chromosome, start, end)


def _coordinates(chromosome, start, end):
    if not isinstance(chromosome, str) or not chromosome:
        raise ValueError("A nonempty chromosome string is required")
    if (
        any(
            isinstance(x, bool) or not isinstance(x, Real) or not isfinite(x) or int(x) != x
            for x in (start, end)
        )
        or start < 1
        or end < start
    ):
        raise ValueError("Coordinates must be finite positive 1-based inclusive integers")


def _model(model, identifier, children=()):
    _coordinates(model.chromosome, model.start, model.end)
    if (
        not isinstance(identifier, str)
        or not identifier
        or model.strand not in {"+", "-", ".", "?"}
    ):
        raise ValueError("A nonempty model ID and valid strand (+, -, . or ?) are required")
    for child in children:
        if child.chromosome != model.chromosome:
            raise ValueError("Child chromosome disagrees with parent model")
        if child.start < model.start or child.end > model.end:
            raise ValueError("Child coordinates extend outside parent model")


@dataclass(frozen=True)
class ExonModel:
    exon_id: str
    chromosome: str
    start: int
    end: int
    strand: str = "."

    def __post_init__(self):
        _model(self, self.exon_id)


@dataclass(frozen=True)
class TranscriptModel:
    transcript_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    exons: tuple[ExonModel, ...] = ()
    tags: tuple[str, ...] = ()
    name: str | None = None
    biotype: str | None = None

    def __post_init__(self):
        _model(self, self.transcript_id, self.exons)


@dataclass(frozen=True)
class GeneModel:
    gene_id: str
    name: str
    chromosome: str
    start: int
    end: int
    strand: str
    transcripts: tuple[TranscriptModel, ...] = ()
    exons: tuple[ExonModel, ...] = ()
    biotype: str | None = None

    def __post_init__(self):
        _model(self, self.gene_id, (*self.transcripts, *self.exons))
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("A nonempty gene name is required")


@dataclass(frozen=True)
class RegionAnnotation:
    chromosome: str
    start: int
    end: int
    genome_build: str
    genes: tuple[GeneModel, ...] = ()
    provenance: dict[str, Any] = field(default_factory=dict)
    chromosome_policy: str = "chr"

    def __post_init__(self):
        validate_region(self.chromosome, self.start, self.end, self.genome_build)
        chrom = normalize_chromosome(self.chromosome, self.chromosome_policy)
        for gene in self.genes:
            if normalize_chromosome(gene.chromosome, self.chromosome_policy) != chrom:
                raise ValueError("Gene chromosome disagrees with annotation region")
            if gene.end < self.start or gene.start > self.end:
                raise ValueError("RegionAnnotation contains a nonoverlapping gene")


def union_exons(exons) -> tuple[tuple[int, int], ...]:
    """Merge overlapping or adjacent inclusive intervals, without mutating models."""
    merged = []
    exons = tuple(exons)
    for exon in exons:
        _coordinates(exon.chromosome, exon.start, exon.end)
    for start, end in sorted((int(e.start), int(e.end)) for e in exons):
        if start < 1 or end < start:
            raise ValueError("Invalid exon coordinates")
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return tuple(merged)


def gene_exons(gene):
    return gene.exons or tuple(e for t in gene.transcripts for e in t.exons)


def select_transcript(gene: GeneModel):
    """Prefer explicit MANE/Ensembl canonical tags; otherwise name the fallback."""
    if not gene.transcripts:
        return None, "no_transcript; collapsed gene exons"

    def rank(t):
        priority = 0 if "MANE_Select" in t.tags else 1 if "Ensembl_canonical" in t.tags else 2
        length = sum(b - a + 1 for a, b in union_exons(t.exons))
        return priority, -length, t.transcript_id

    transcript = min(gene.transcripts, key=rank)
    priority = rank(transcript)[0]
    rule = (
        "annotated MANE_Select",
        "annotated Ensembl_canonical",
        "representative fallback: longest union-exon length, then transcript ID",
    )[priority]
    return transcript, rule


def matches_gene(gene, target):
    if target is None:
        return False
    return target in (gene.gene_id, gene.name) or gene.gene_id.split(".")[0] == target.split(".")[0]
