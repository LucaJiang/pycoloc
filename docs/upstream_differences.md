# Pinned upstream differences and correctness checkpoint

Reference: coloc 6.0.3, commit 8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45,
archive SHA256 f6580ab471e1bd9f3acc45e23a73bb2ee16fe3c0a45e01cd8cc008854661b479.
These findings concern that exact source. No R outputs were edited. Version
0.1.1 makes mathematically corrected fields the default (result schema 2).

## Confirmed implementation errors

**COLOC-WEIGHTED-SNP-001 — claudia.R coloc.abf; susie.R coloc.bf_bf/coloc.susie.**
The upstream main SNP.PP.H4 omits weights despite weighted H4 evidence. With
BF1=BF2=[1,1], p1=p2=0.0001, p12=0.000005, w1=[9,1], w2=[1,1],
q12=[0.000009,0.000001]. R and Python 0.1.0.dev0 main values are [0.5,0.5];
the old Python .weighted field and corrected main are [0.9,0.1]. With both
weights [9,1], q12=[0.0000162,0.0000002], sum=0.0000164, and the corrected
conditional is [81/82,1/82]. No fixed-sum q12 renormalization is allowed.
Independent Decimal configuration enumeration tests
`test_analytic_weighted_main_all_interfaces` cover ABF/BF/SuSiE and both engines;
`test_identity_symmetry_every_pair_explicit_null` checks every pair.
.weighted is now an equal alias. Existing R micro and real outputs remain as
legacy_R columns in new-run audit TSVs, alongside independent oracle values.

**COLOC-WEIGHT-ORDER-001 — claudia.R merge/membership; susie.R column-count shortcut.**
For SNPs c,a,b with nonuniform ID-bound weights, upstream ABF sorts the merged
SNPs but subsets weights in original input order. BF pairing omits weight
reindexing when complete-overlap column counts match. Python's existing ID
binding is retained. Executed historical four-SNP probes recorded max H0–H4
errors 0.00008614739250312198 (ABF) and 0.14369617072589103 (BF). Original R and
Python values and exact fixtures are in the old A_micro.json and retained
pinned-R audit records. Current tests independently enumerate the
expected configurations, permute either/both inputs with bound weights, scale
weights and exchange traits, with equal column counts and partial intersection:
`test_abf_independent_permutation_and_trait_exchange` and
`test_identity_symmetry_every_pair_explicit_null`. The latter includes explicit
null offsets and a biological SNP literally named null.

**COLOC-EFFECT-PRIOR-001 — claudia.R coloc.abf/process.dataset.**
The R wrappers accept but drop overrides before approx.bf.estimates. At beta=0.6,
V=0.04, cc SD=0.3 or 0.7, the correct logABF is the log normal density ratio
N(beta;0,V+SD²)/N(beta;0,V); R uses SD=0.2 instead. Quantitative W=(SD*sdY)²,
not SD*sdY or SD itself. Python already forwarded these overrides and continues
to do so. The executed old four-SNP cc override (cc=0.4) has max hypothesis
error 0.027988814116981064 against R. `test_effect_prior_normal_integral_public`
independently checks two nondefault SDs for both traits and both output ABFs.
The pinned source documentation calls effect_priors variances, while its code
squares them; public Python docs now state SD units explicitly.

## Numerical defects and stable equivalents

**COLOC-H3-001 — claudia.R logdiff of product and diagonal masses.**
For logBF1=logBF2=[1000,0], R returns an off-diagonal log mass of −∞ because both
inputs to subtraction round to 2000. The exact unprior-weighted value is
1000+log(2)=1000.6931471805599; logU3 additionally includes log(p1)+log(p2).
Python's existing positive O(m) summation is unchanged. Tests
`test_h3_extreme_coincident_lead_does_not_cancel` and
`test_zero_extreme_and_underflow_public` check both internal log evidence and
final hypotheses with independent enumeration. Old R output is retained in
the retained private H3 reference audit. ABF shrinkage likewise
retains stable log1p/logaddexp evaluation rather than log(1-r) after r rounds to 1.

**COLOC-ZERO-SUPPORT-001 — all-negative-infinity logsum.**
For logBF1=[−∞,0], logBF2=[0,−∞], R logsum of shared terms becomes NaN. The
public BF wrapper subsequently returns H0=1 and H1=H2=H3=H4=0, plus NaN
conditional SNP PP (actual revision probe: the retained private official probes).
That fallback loses positive H1/H2/H3 evidence. Python retains H4=0, H3>0 and NaN conditional
SNP PP. Version 0.1.1 additionally reports logU4=−∞ and
h4_status=zero_shared_evidence. Disjoint prior support has the same meaning.
For logBF1=[−1000,−1001], logBF2=[0,0], finite logU4 and underflowed PP.H4=0
instead give the valid conditional softmax([0,−1]) and
h4_status=posterior_underflow. This is diagnostic refinement, not a change to
Python's existing log-space conditional formula. See
`test_zero_extreme_and_underflow_public`, `test_abf_h4_underflow_is_defined`,
`test_each_pair_has_its_own_zero_or_underflow_state` and
`test_all_negative_infinity_does_not_contaminate_hypotheses`.

## API-domain and representation choices

**COLOC-ZERO-WEIGHTS-001 — nonnegative biological weights.**
Pinned ABF permits partial zeros; BF/SuSiE reject them through <=0 validation.
Python previously copied that distinction. All three interfaces now accept
[1,0] with positive per-side total. With unit BFs, w1=[1,0], w2=[0,1], U3 is
positive, U4=0 and conditional SNP PP is NaN. This is legitimate input-domain
support, not an allegation of an upstream fitting algorithm error. Empty
post-intersection support raises ValueError; no uniform fallback or epsilon.
Tests: analytic weighted fixture above, `test_empty_intersection_weight_support_errors`,
and the retained ABF `test_no_weight_mass_on_intersection`.

**COLOC-PVALUE-PRIOR-001 — explicit effect_priors on P/MAF/N.**
Pinned approx.bf.p exposes fixed prior SDs 0.15 (unit-variance quant) and 0.2
(logOR cc), with no override argument. R and old Python silently used these
scales when the wrapper received a mapping. Python 0.1.1 rejects that combination
with ValueError, for either trait and either input side; the default P branch
and beta/variance precedence remain unchanged. This is an explicit API boundary
and a fix to Python's silently ignored argument, not a new ABF model.
Tests: `test_explicit_effect_prior_pvalue_combination_rejected` and retained
`test_pvalue_formula`/`test_abf_formula_and_beta_precedence`.

**Null, positions, overlap and structural states (unchanged).** Null is declared
by null_index=m and excluded from biological IDs/weights. A biological SNP named
null remains biological; R uses that name as a sentinel. ABF matches IDs and
rejects conflicting positions instead of merging on all identically named
columns. No-CS, no-overlap and trimming remain distinct states. Overlap is the
BF-only posterior fraction on common biological support conditional on an
association, computed before intersection, without coloc variant weights.
Fine-mapping prior weights are never reapplied to BF rows, and alpha is never
used in place of BF. These are documented API/representation choices, not blanket
upstream-error exemptions. Scalar positive base priors remain required.

## Executed evidence locations

Historical and correctness audit records are retained in the private development evidence. Public analytical regressions below reproduce the mathematical counterexamples.

## 0.2.0 visualization extension

`plot_coloc_locus`, local annotation models/indexing, separate per-trait LD colors
and the vector genomic LD triangle are new Python plotting interfaces, not an
upstream numerical discrepancy. Classification: API extension. No statistical
formula, main-field correction, result schema, tolerance or documented 0.1.1
upstream difference changes. Model/coordinate/index/LD validation errors are
explicit visualization input boundaries; no new claims of upstream defects.

## 0.2.1 equivalent implementation changes

Classification: implementation/allocation optimization, with no new numerical
difference ID. Reusing a shared log-sum and building paired result columns once
does not change shared q12, stable H3 or zero/underflow semantics. Existing
correctness difference IDs and frozen R bounds remain applicable. The final
installed-wheel comparison includes both engines, independent oracles and exact
before/after result comparisons; timings are separate from correctness evidence.
