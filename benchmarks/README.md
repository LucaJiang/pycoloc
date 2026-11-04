# Full-call comparisons with R coloc

These manual scripts compare ABF, all supplied component BF rows, or existing
SuSiE fit pairing. They use the installed pycoloc public API and pinned coloc
6.0.3. They do not fit SuSiE; fit-to-fit numerical differences are assessed
separately in the [agreement report](../docs/r_agreement.md).

Run from the checkout after installing pycoloc. R additionally needs jsonlite,
coloc 6.0.3 and susieR 0.16.6 from the pinned sources in
[reference mapping](../docs/reference_mapping.md). Set `R_LIBS_USER` to that
library. Use an optimized R BLAS, verify its actual thread count, and record the
R/package source identity and CPU/OS when comparing performance.

## Inputs and numerical checks

Export one complete bundled toy pair into a **new** directory:

```sh
python benchmarks/coloc_case.py --export-fixture bench-output/D1_D2 --pair D1 D2
Rscript --vanilla benchmarks/coloc_case.R bench-output/D1_D2/case.json abf validation bench-output/R-abf.json
python benchmarks/coloc_case.py --case bench-output/D1_D2/case.json --stage abf --mode validation --reference bench-output/R-abf.json --output bench-output/Python-abf.json
```

Use `--pair D3 D4` for the other complete 500-SNP pair. The exporter copies the
stored float64 input/BF values bit for bit into little-endian files, records
shapes and SHA256 hashes, and retains SNP order and original component IDs in
`case.json`. No LD is loaded for colocalisation. Both implementations load only
the data for the requested stage.

Change `abf` to `bf` or `susie` in both commands for the other stages. `bf` uses
**all supplied BF rows**; `susie` calls the public fit consumer, which selects
credible-set component IDs. Keep these workloads distinct. The ABF prior is
p12=1e-5; both pairing priors are p12=5e-6, with p1=p2=1e-4. Parameters are
explicit in `case.json` and never inferred from the timing outcome.

The comparison checks H0–H4 with atol=1e-10, rtol=0 and conditional SNP
probabilities/logBFs with atol=1e-10, rtol=1e-10. IDs, order and structural states
compare exactly. No-CS/no-overlap remain unavailable states; failed fields are
retained. Intentional corrected-weight/numerical cases have independent oracles
in the existing tests and `examples/check_example.py`.

## Timing and peak process RSS

Run timing without other heavy work. Set one numerical thread in both processes:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p bench-output/clock
cp benchmarks/monotonic_clock.c bench-output/clock/
(cd bench-output/clock && R CMD SHLIB monotonic_clock.c)
Rscript --vanilla benchmarks/coloc_case.R bench-output/D1_D2/case.json abf timing bench-output/R-time.json 7 10 bench-output/clock/monotonic_clock.so
python benchmarks/coloc_case.py --case bench-output/D1_D2/case.json --stage abf --mode timing --repeats 7 --batch 10 --output bench-output/Python-time.json
```

Each backend performs one warmup and seven retained batches of ten full calls.
Monotonic wall clocks include public validation, preprocessing, computation and
standard result construction, excluding imports, input loading and numerical
export. Only one returned object is retained at a time. R's public-call status
printing is captured in memory and exported after timing and peak RSS capture;
warnings remain visible. Alternate
backend order by predeclared case index in a multi-case comparison. Batch start
and end epoch times support checking external interference.

`--engine numba` is optional. Its first warmup includes lazy import/compilation;
record that separately from its warm retained calls. Default NumPy is the primary
comparison. Repeated calls always recompute; there is no fitted-output cache.

Measure memory in separate, fresh processes, once per case/backend:

```sh
Rscript --vanilla benchmarks/coloc_case.R bench-output/D1_D2/case.json abf memory bench-output/R-memory.json
python benchmarks/coloc_case.py --case bench-output/D1_D2/case.json --stage abf --mode memory --output bench-output/Python-memory.json
```

On Linux, **peak process RSS** includes interpreter/imports, input loading,
allocations during the full call and its returned object. It is captured before
serializing large numerical results. Python uses `ru_maxrss/1024`; R uses Linux
`VmHWM/1024`, both in MiB. This is an absolute process high-water mark, not an
estimate of inference-only memory. Timing runs intentionally leave memory
metrics null. Errors have explicit messages rather than invented zero timings.

Use median R time / median Python time as per-case speedup; peak R / peak Python
as the separate memory ratio. Report denominators, wins/ties/regressions, IQR and
all unavailable or nonconverged cases. Freeze case selection before results.
Full real-data measurements are manual and outside default CI.

## Adapters and checks

`load_case(path, stage)` loads this exchange schema. `run_case(inputs, stage)`
returns an unmodified `ColocResult`, so an external driver can supply already
loaded, audited real inputs. `compare_result(result, raw_R_result)` checks ABF or
pairing; pass raw R's one-based original component IDs and explicit unavailable
status. The same comparison can assess native-fit → pycoloc against R-fit →
coloc without treating PIP agreement as posterior agreement.

```sh
python -m pytest benchmarks/test_coloc_case.py -q
```
