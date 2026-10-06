# ColocResult schema 2 (pycoloc 0.1.1)

`metadata.result_schema_version` is the string `2`. Summary shape and all
existing H0–H4, prior, SNP and component identifiers retain their definitions.
`SNP.PP.H4` for ABF, and the pair-specific column named by
`summary.snp_pp_column` for BF/SuSiE, is P(shared causal SNP=j | H4, data).
It uses the exact q12 from hypothesis evidence. Each `.weighted` column is an
identical alias; there is no implicit legacy/unweighted output mode.

Computed summary rows add these diagnostics:

| Field | Definition |
|---|---|
| logU3 | Natural log of off-diagonal configuration evidence; −∞ for zero support |
| logU4 | Natural log of shared configuration evidence; −∞ only for zero support |
| h4_status | defined, zero_shared_evidence, or posterior_underflow |

True U4=0 returns PP.H4=0 and all-NaN conditional SNP posteriors. Finite logU4
with PP.H4 underflowing to zero retains a conditional distribution summing to 1.
Partially-zero weights give exact zeros in the conditional when defined.
An input side with zero total weight after intersection raises ValueError.

`result.status` continues to describe structural eligibility: ok,
partially_trimmed, insufficient_overlap, no_overlap, or no_cs. Rows retained only
as insufficient-overlap diagnostics do not claim log evidence or h4_status.
These are not five-hypothesis posteriors. No-CS and no-overlap remain empty
results, distinct from mathematical U4=0.

Candidate tables should sort the main column for the selected summary row.
Multiple pairs must retain their original idx1/idx2; never aggregate them into
a region-level H4. A shared joint posterior, if needed, is PP.H4 times the main
conditional SNP posterior, not the conditional alone.


## Result container and rows

`ColocResult` exposes summary/results pandas DataFrames, priors, method, status,
weights and metadata. Attribute access and result['summary'] are equivalent.
as_dict() is shallow; nsnps is the common biological-SNP count (zero for empty
summary). metadata.result_schema_version is "2". weights records common-SNP
indexed weights. metadata retains convergence, original support and credible sets.

ABF has one computed summary row. BF/SuSiE pair rows include nsnps, PP.H0.abf
through PP.H4.abf, original idx1/idx2, corresponding idx1_r/idx2_r, hit1/hit2,
overlap1/overlap2, status and snp_pp_column. Row order can change when input
components/credible sets are returned in a different order; use component IDs.
The pair's snp_pp_column names its main conditional posterior in results; never
guess the column from a retained-set ordinal. With insufficient overlap,
diagnostic rows are not valid five-hypothesis distributions.

Result tables preserve SNP identifiers and relevant association/BF fields. No
scalar H4 is synthesized from multiple signal pairs. Nonconverged inputs emit
a warning and remain visible in metadata; no refit or hidden R execution occurs.
