# Bayesian colocalisation

<p class="intro">Compare association signals with Bayesian colocalisation in Python. Start with summary statistics or existing SuSiE fits, then inspect the hypotheses, shared-variant probabilities, and model assumptions.</p>

<div class="cards">
<a class="card" href="install.html"><span>01 / Start</span><strong>Install pycoloc</strong><small>Pure Python with NumPy and SciPy. R is not a runtime requirement.</small></a>
<a class="card" href="example.html"><span>02 / Verify</span><strong>Run the 500-SNP check</strong><small>Frozen public teaching data, official-R expectations, and useful failure messages.</small></a>
<a class="card" href="r_agreement.html"><span>03 / Understand</span><strong>Read the R agreement report</strong><small>Executed comparisons, absolute errors, and intentional corrected-field differences.</small></a>
</div>

## Quick start

Follow [installation](install.md) to clone the repository and install its dependencies, then run from the checkout:

```sh
python examples/check_example.py --output-dir example-output
```

Once dependencies are installed, the check needs no network, R, private study files, or prusie. It calls real public APIs on the bundled upstream **synthetic teaching fixture**. It checks file integrity and numerical results, writes a machine-readable record, and exits nonzero on a substantive mismatch.

## Choose your input

| Input you have | Public function | Interpretation |
| --- | --- | --- |
| Two sets of beta and variance, or P values with required metadata | `coloc_abf` | Single causal variant per trait; no LD required |
| Full component log Bayes-factor matrices | `coloc_bf_bf` | A result for each component pair |
| Existing SuSiE results with credible sets | `coloc_susie` | Pairs the original credible-set component rows |

Both coloc prior-weight vectors follow explicit SNP identities. An existing fit is consumed directly; fine-mapping is never run implicitly. Read [inputs and alignment](inputs.md) before preparing an analysis.

## Read the result in context

`PP.H0.abf` through `PP.H4.abf` describe the five colocalisation hypotheses for one comparison. `SNP.PP.H4` is conditional on H4: it is not the unconditional probability of sharing a causal SNP. A multiple-signal result has one row per component pair, without an invented region-wide H4.

<div class="callout"><p><strong>Agreement has a defined scope.</strong> Unweighted paths with explicitly recorded priors are checked against pinned official R. Corrected weighted SNP posteriors, stable H3, and several input semantics deliberately differ where the reference is inconsistent. The report preserves those differences and checks independent mathematical expectations.</p></div>

The [result schema](result_schema.md), [methods and priors](methods.md), and [upstream differences](upstream_differences.md) explain how to interpret computed, missing, and unavailable results. [Locus plotting](locus_plot.md) adds association tracks and local annotation without changing inference.

## Reference

The algorithm reference is official coloc 6.0.3; the example's SuSiE fits use
susieR 0.16.6. See [compatibility](compatibility.md) for supported scope and
[source mapping](reference_mapping.md) for pinned identities.
Maintained by Wenxin Jiang. Source: [GitHub](https://github.com/LucaJiang/pycoloc).
