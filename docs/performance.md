# Performance and resource use

The NumPy hypothesis implementation is O(m) and avoids an m×m configuration
matrix. Optional Numba accelerates the numerical hypothesis kernel; report
first compilation separately from warmed calls. R and fine-mapping are never
invoked implicitly. Independent-region process concurrency is controlled by the
caller, who must budget total library threads and memory.

LD heatmap thinning affects only visualization. It reindexes both matrix axes
by one explicit SNP permutation, records omissions, and preserves association
tracks and inference inputs. Stored SNP order need not be genomic order.
Whole-genome GTF indexing can take substantial time and disk space; rendering
prepared local annotation avoids rebuilding the index. Cache identities retain
source, release, build, indexer/schema and content hashes.

Prusie comparisons depend on actual R BLAS/LAPACK and LD redundancy; no such
timings establish an independent pycoloc speed claim. The package supplies no
whole-host utilization guarantee or hidden worker pool. See the separate
prusie performance guide when fine-mapping is part of your analysis.
