# Development check notes for fastFGEE 0.3.0.9006

This is a development candidate reconstructed from the validated 0.3.0.9004
source tree. It must not be submitted to CRAN until the commands in
`validation/run_R_validation.sh` have completed successfully in an R-enabled
environment.

## Package-controlled issues addressed

* The archived `irregulAR1` dependency was removed. The package contains an
  independent exact tridiagonal precision implementation following Allevius
  (2018).
* The `sanic` dependency was removed. Symmetric positive-definite operations use
  registered Rcpp/LAPACK routines with a base-R fallback.
* Rcpp is a direct `Imports` and `LinkingTo` dependency; `sourceCpp` is not
  imported or used.
* All Rcpp entry points are registered.
* `refund (>= 0.1-40)` is declared.
* The public estimator surface is one-step only.
* The stray `Rplots.pdf` artifact is excluded and checked by the bundle audit.

## Check status in the construction environment

R and Rscript were not installed in the environment that assembled this source
candidate. Consequently, no claim of a successful `R CMD check` is made here.
The bundle records the static audit and independent numerical checks and
contains a fail-fast script for build, install, tests, and `--as-cran` check.
This file should be replaced with the actual final check environments and
results after that script passes.
