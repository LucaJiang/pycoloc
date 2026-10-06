# Methods and parameters

The ABF uses a normal effect prior with standard deviation 0.2 for log odds
ratios, or 0.15 × sdY for quantitative beta estimates. For P-value input the
quantitative prior SD is 0.15, using variance 1 / (2 N f (1−f)); case-control
variance divides additionally by s(1−s). With effect variance W and sampling
variance V, r = W/(W+V), and logBF = ½[log(V/(V+W)) + r z²]. All arrays use
float64 and natural log BF. Input P values are two-sided.

`effect_priors={'quant': ..., 'cc': ...}` specifies **standard deviations**,
not variances. On beta/varbeta inputs quantitative SD is multiplied by sdY
before squaring; case-control SD is on the log(OR) scale. Nondefault values
are forwarded to the actual ABF calculation. The active P/MAF/N branch retains
fixed SDs 0.15 (unit-variance quantitative phenotype) and 0.2 (logOR), and
rejects an explicit effect_priors mapping. Supplying both input branches still
selects beta/varbeta and applies the mapping.

When sdY is unavailable, regress 2 N f (1−f) against 1/varbeta through the origin
and take the square root of its slope. Record the estimate in
`metadata['dataset1'/'dataset2']`; the software warns about estimation. N is
scalar. Frequency validation matches pinned coloc's open interval (0,1);
the variance expression is invariant to f versus 1−f, but frequencies must not
be used to infer allele direction without effect-allele evidence.

For Q common SNPs, weights define q1[j]=Q p1 w1[j]/Σw1 and analogously q2.
Omitting a side's weights means q=p for that side. Shared q12[j] is
q1[j]q2[j]p12/(p1p2), **without** fixing its sum to Q p12. H1 and H2 sum
single-trait contributions; H3 sums q1[i]q2[j] BF1[i]BF2[j] over i≠j; H4 sums
q12[j] BF1[j]BF2[j]. H0 has unnormalised mass 1. Normalize the five masses.

H3 uses O(Q) log-prefix/suffix reductions (NumPy) or an O(Q) sequential
logaddexp kernel (Numba), avoiding subtraction of nearly equal masses and any
Q×Q allocation. `engine='numba'` compiles on first call, without fastmath or a
disk compilation cache. The default is the readable NumPy implementation.

Base priors must be finite scalars in (0,1]. Boundary zero priors are rejected
explicitly: pinned weighted formulas divide by p1p2 and do not define a reliable
zero-prior convention. Vector base priors are unsupported; use SNP weights.
ABF adjusts p1 and p2 using the original trait SNP counts and p12 after
intersection: when n p≥1, replace p by 1/(n+1), with a warning. BF pairing does
not apply that adjustment to hypothesis priors. Weight normalization occurs on
the common support after these steps.

Weights may be keyed mappings/Series or arrays. Arrays first bind to the
input's explicit SNP IDs. Keys must match each complete biological support.
All three interfaces accept partial zeros and reject all-zero common support.
Infinite, missing, negative or mismatched weights always raise. Zero weights
remove configuration support exactly, without an epsilon or uniform fallback.
The main SNP.PP.H4 fields normalize q12[j] BF1[j] BF2[j] in log space;
`.weighted` fields are equal transitional aliases. No extra q12 normalization
is applied. This is P(shared SNP=j | H4, data); a joint posterior would multiply
it by PP.H4. Disjoint shared support gives logU4=−∞, PP.H4=0, conditional NaN,
and h4_status='zero_shared_evidence'. Finite logU4 with PP.H4 rounded to zero
gives h4_status='posterior_underflow' and a valid conditional posterior. Ordinary
computed rows have h4_status='defined'. Both logU3 and logU4 are exposed for
diagnostics. These statuses are separate from structural result.status.

SuSiE consumption selects original `sets.cs_index` **rows** of lbf_variable and
retains every biological column, including SNPs outside CS membership. Rows
pair with trait 1 varying fastest. Credible set indices and all public component
indices are zero based; idx1_r/idx2_r explicitly add one for R comparisons.
An optional null must be the final column and be declared by null_index=m.
Biological IDs and coloc weights exclude it. Each row is shifted by its null
logBF. Overlap is the posterior fraction on common SNPs conditional on a
biological association, evaluated against each full input support. It ignores
coloc weights just as the pinned source does. It is not a SNP-count overlap.

`no_cs`, `no_overlap`, `insufficient_overlap`, `partially_trimmed` and `ok` are
distinct. Insufficient overlap for every pair preserves rows with H1–H4 missing
and the upstream minimum-null H0 diagnostic; this is not a valid five-hypothesis
posterior. Unconverged fits emit a warning and retain explicit convergence
metadata. Computation is still possible, matching upstream consumption.

`plot_coloc` draws full supplied regional marginal tracks. Original valid P
values take priority over Wald P. Natural logP or stable normal log tails recover
underflowed zeros when possible; otherwise clipping is labelled. No conditional
association or LD coloring is inferred. Coordinates and chromosome are mandatory.
Titles derive H0–H4 from one exact summary row. Multiple pairs require explicit
original component IDs. The figure and axes remain editable; save as SVG, PDF
or PNG without a GUI.
An outlined diamond marks the shared candidate chosen from the corrected main
SNP posterior for that exact row. It is absent when shared support is undefined.
Component leads continue to identify the maximum component BF separately.

See [performance and resource use](performance.md) for JIT, annotation and plotting costs.
