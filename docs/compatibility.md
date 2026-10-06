# Compatibility scope

Supported: beta/varbeta and P-value/MAF/N ABFs; cc and quant traits; known and
estimated sdY; scalar base priors; both coloc weight vectors in every interface;
all-pair BF computation; full posterior overlap and trimming; explicit null;
existing SuSiE results with nonconsecutive CS indices; marginal regional plots.

The algorithm reference is pinned coloc 6.0.3, official GitHub commit
8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45, source archive SHA256
f6580ab471e1bd9f3acc45e23a73bb2ee16fe3c0a45e01cd8cc008854661b479.
The archive hash in this document is checked against the project source lock
during validation; consult the source mapping for provenance.

The exported SusieResult contract is a mapping or object with `as_dict()`.
Required consumption fields are `variant_ids`, `lbf_variable`, `sets` and
`converged`; null_index is explicit when present. Additional arrays and metadata
are allowed. pycoloc never imports prusie during ABF or fit consumption.

Not implemented: raw-dataset automatic SuSiE fitting, runsusie retries,
back-calculation of BF from alpha, finemap APIs, conditioning/masking,
coloc.signals, sensitivity-analysis functions, vector base-prior legacy paths,
or RDS import. The newer `plot_coloc_locus` interface supports explicitly supplied LD colors; `plot_coloc` remains a marginal-track plot. Unknown keyword arguments fail rather than silently
being ignored. Runtime calls to R are absent. This is not a generic backend or
workflow framework.

Intentional changes are detailed in upstream_differences.md: keyed weights,
explicit null support, stricter finite-value validation, position mismatch
errors, stable H3 and log shrinkage calculations, uniform result objects,
explicit nondefault effect-prior support on beta/varbeta, explicit rejection
on P/MAF/N, and weighted main SNP posterior fields with equal .weighted aliases.
These must not be described as exact parity on upstream-bug counterexamples.

Version 0.1.1 uses ColocResult schema 2 (metadata.result_schema_version).
Partial-zero weights are now legal for every coloc interface; they denote exact
zero configuration support. The overlap diagnostic keeps its previous BF-only
definition and is unaffected by coloc weights. See [result schema](result_schema.md).
