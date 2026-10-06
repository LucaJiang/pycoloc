"""Offline annotation models and optional local GTF/GFF3 indexing."""

from .local import annotation_from_dict, annotation_to_dict, prepare_region_annotation
from .models import ExonModel, GeneModel, RegionAnnotation, TranscriptModel, union_exons

__all__ = [
    "ExonModel",
    "GeneModel",
    "RegionAnnotation",
    "TranscriptModel",
    "prepare_region_annotation",
    "annotation_from_dict",
    "annotation_to_dict",
    "union_exons",
]
