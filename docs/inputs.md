# Inputs and statistical interpretation

For ABF, each dataset supplies unique SNP IDs and trait type `cc` or `quant`.
Use signed beta (log OR for case-control) and positive varbeta=SE². Quantitative
traits need known sdY or MAF/N for the documented estimator. The alternative
P-value branch needs P, MAF, N and case fraction s for case-control. If both
branches are present, beta/varbeta takes precedence. Explicit effect_priors
are prior SDs and apply only to that branch. ABF requires no LD.

Weights can be SNP-keyed mappings/Series or vectors bound to the supplied SNP
order. Each side must have positive total after intersection; partial zeros are
allowed. For Q common SNPs, q1=Q*p1*w1/sum(w1), q2=Q*p2*w2/sum(w2) and
q12=q1*q2*p12/(p1*p2). q12 is not separately normalized. Main SNP.PP.H4 includes
both weights; the `.weighted` columns are equal aliases. True zero shared
evidence yields undefined conditional SNP posteriors; H4 underflow can still
have a valid conditional distribution. See [result schema](result_schema.md).

`coloc_bf_bf` consumes full natural-log BF rows. `coloc_susie` selects the full
lbf_variable rows named by original sets.cs_index, preserving all biological SNP
columns, not only CS members. Marginal PIPs are not substitutes for component
BFs. Fine-mapping prior_weights and coloc prior_weights1/2 act at different stages;
do not apply one again automatically. Null columns must be declared explicitly.

H0–H4 are comparison-specific hypotheses. A conditional SNP posterior is
P(shared SNP | H4,data), not the joint H4/SNP probability. Multiple component
pairs do not define one aggregate region-level H4. no_cs, no_overlap,
insufficient_overlap and nonconvergence are distinct states; unavailable
probabilities must not be replaced by zeros.

For fine-mapping with prusie, use signed r, exact SNP/allele order and documented
LD population/build/sample provenance. Plotting accepts r or r² only when kind
is explicit; r² is a display quantity. Heatmap row/column selection is coordinated
by SNP IDs and positions and never changes the association inputs or inference.
Local annotation, coordinate conventions and cache metadata are documented in
[locus plots](locus_plot.md). Only the synthetic toy example is bundled; real genotype, study summary
statistics and annotation databases must be supplied by the caller.
