# Differences from pinned coloc

These differences concern coloc 6.0.3 at commit
`8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45`. pycoloc preserves corrected
mathematical behavior rather than reproducing inconsistent upstream fields.
The [executed report](r_agreement.md) keeps raw R results separate from
independent mathematical expectations; [source mapping](reference_mapping.md)
identifies the relevant functions.

## Weighted SNP probabilities and alignment

The reference includes variant weights in H4 evidence but omits them from its
main `SNP.PP.H4` conditional probabilities. With unit Bayes factors and weights
[9,1] and [1,1], R returns [0.5,0.5]; the shared-prior posterior is [0.9,0.1].
pycoloc's main field uses the same q12 as its H4 evidence. `.weighted` is an
identical alias. With both weight vectors [9,1], the correct conditional is
[81/82,1/82]. Shared q12 is never renormalized to a fixed total.

The reference ABF merge can sort SNPs without reordering weights, and its BF
path can skip reindexing when overlap column counts match. pycoloc binds each
weight to its input SNP ID before sorting or intersection. Regression tests
independently enumerate configurations and permute either trait, including
partial intersection and explicit null columns. See
`test_analytic_weighted_main_all_interfaces`,
`test_abf_independent_permutation_and_trait_exchange` and
`test_identity_symmetry_every_pair_explicit_null` in the public tests.

## Effect-prior scales

On the beta/varbeta route, pycoloc honors explicit `effect_priors` while the
pinned R wrapper drops those overrides. Values are standard deviations:
quantitative W=(SD × sdY)², and case-control W=SD² on the log-OR scale. The
independent `test_effect_prior_normal_integral_public` checks the corresponding
normal-density ratio for both trait types and multiple SDs.

The P/MAF/N route uses fixed SDs 0.15 for unit-variance quantitative phenotypes
and 0.2 for log OR. It rejects explicit `effect_priors` instead of silently
ignoring them. Supplying both routes still selects beta/varbeta. See
`test_explicit_effect_prior_pvalue_combination_rejected` and [methods](methods.md).

## Stable H3 and shared support

For logBF1=logBF2=[1000,0], subtracting total and diagonal masses in the
reference loses the off-diagonal mass. Its unweighted log value is actually
1000 + log(2). pycoloc evaluates the positive off-diagonal sum in O(m) time,
without an m×m configuration matrix or cancellation. Stable ABF shrinkage also
avoids evaluating log(1−r) after r rounds to one.

For logBF1=[−∞,0] and logBF2=[0,−∞], shared evidence is zero but distinct-causal
H3 evidence remains. pycoloc preserves it; the pinned R wrapper falls back to
H0=1 after an invalid shared log sum. Disjoint prior support has the same
zero-shared-evidence meaning. Tests include
`test_h3_extreme_coincident_lead_does_not_cancel` and
`test_all_negative_infinity_does_not_contaminate_hypotheses`.

True zero shared evidence gives `logU4=−∞`, H4=0 and an undefined conditional
SNP posterior. Finite logU4 with H4 underflowing to zero still has a defined
conditional distribution. `h4_status` distinguishes these cases; see
[result schema](result_schema.md) and `test_abf_h4_underflow_is_defined`.

## Input and representation boundaries

All three interfaces accept partially zero, nonnegative weights with positive
common total; R's BF/SuSiE wrappers reject zeros. Zero weights remove support
exactly, without epsilon or uniform fallback. With unit BFs and disjoint
weights [1,0] and [0,1], H3 remains positive and H4 is zero. All-zero support
after intersection raises an error.

Null columns are explicitly declared with `null_index=m` and excluded from
biological IDs/weights. A biological SNP named `null` remains biological.
ABF matches IDs and rejects conflicting positions. SuSiE pairing uses full BF
rows indexed by original credible-set component IDs, not marginal PIPs or
only set members. Fine-mapping priors are not reapplied to BF rows.

No-CS, no-overlap, insufficient-overlap and nonconvergence retain distinct
states. The overlap measure is BF-based posterior support on common biological
variants, conditional on association, before coloc weighting. These API choices
and the Python plotting interfaces are documented extensions, not claims of
blanket compatibility with every R API.
