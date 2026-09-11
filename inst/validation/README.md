# Validation and performance checks

These development scripts compare the optimized engine with the preserved CRAN
0.1.0 path. They are intended for a workstation or HPC compute node, not for
CRAN examples.

For the selective `0.3.0.9006` reconstruction, start with the fail-fast driver
shipped in the complete development bundle:

```bash
bash validation/run_R_validation.sh /path/to/fastFGEE_0.3.0.9006
```

It regenerates the Rcpp registration files in a disposable copy and requires
an exact diff, then runs `R CMD build`, installation, the focused numerical
checks, the full testthat directory, and `R CMD check --as-cran --no-manual`.
Its output directory is created beside the source tree by default, never inside
the package being validated.

## Prerequisites

Install the development package and its suggested validation dependencies,
including `testthat`, `refund`, `mgcv`, and `SuperGauss`. Rcpp is a direct
package dependency and is installed with fastFGEE.
Use one BLAS/OpenMP thread per independent benchmark process unless the purpose
is explicitly to benchmark threading.

## Recommended order

1. `Rscript inst/validation/run_golden_checks.R`
2. `Rscript inst/validation/stress_working_stats.R`
3. `Rscript inst/validation/benchmark_correlation_operators.R`
4. `Rscript inst/validation/benchmark_working_stats.R`
5. `Rscript inst/validation/benchmark_tuning.R`
6. `Rscript inst/validation/benchmark_fgee_engines.R`
7. `bash inst/validation/run_validation.sh`
8. `bash inst/validation/run_package_checks.sh`

`run_validation.sh` wraps each R process with GNU `/usr/bin/time -v` when
available, so its logs contain process-level peak RSS. CSV outputs are written
to `inst/validation/results/` by default. Override that location with:

```bash
export FGEE_VALIDATION_OUT=/path/to/results
```

The benchmark scripts accept environment variables for larger stress tests.
For example:

```bash
export FGEE_BENCH_N=100
export FGEE_BENCH_NI=100
export FGEE_BENCH_L=50
export FGEE_BENCH_REPS=5
bash inst/validation/run_validation.sh
```

The exact variables used by each script are documented near the top of that
script. Start small, verify numerical equivalence, and then increase dimensions.

## Required acceptance checks

The release candidate should demonstrate:

- regular and irregular AR(1) and exchangeable operators agree with dense
  references within stated tolerances, and regular-grid operators agree with
  the optional SuperGauss reference;
- batched Kronecker results agree with dense calculations;
- one-pass `W_i,d_i` agree with the legacy separate constructors;
- exact-Gaussian `d_i + W_i beta0` agrees with the legacy third pass;
- staged Gaussian sufficient-statistic fastK reproduces the current objective;
- analytic fastK and qREML gradients agree with finite differences;
- compact sandwich and wild-cluster results agree with legacy list-based code;
- optimized and legacy one-step fits agree at fixed seeds and grids;
- peak memory and elapsed time are reported separately for working-statistics
  construction, tuning, and the entire fit.

The optimized package should not be promoted solely on microbenchmarks. Run
fresh-seed simulations and at least one application-sized dataset before a
release.
