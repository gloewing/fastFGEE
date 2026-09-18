## Submission

This release does two things: it fixes the test failures reported on the CRAN
check page for 0.2.0, and it adds a working-statistics optimisation plus a
revised vignette.

## The failure being fixed

Six check flavors reported:

    Error in `.source_root()`: Could not locate the fastFGEE source tree
    for source-level tests.
    test-dependency-and-registration.R:22, :53, :95
    [ FAIL 3 | WARN 2 | SKIP 6 | PASS 1969 ]

Three tests in `tests/testthat/test-dependency-and-registration.R` are
source-level audits: they read `DESCRIPTION`, `R/*.R` and `src/*.cpp` directly
to confirm that no archived dependency is referenced, that every Rcpp attribute
export is registered, and that each core wrapper has a single active definition.
Those files do not exist when the suite runs against an installed package, which
is the case on the affected flavors, so the helper that locates the source tree
raised an error.

The helper now returns `NULL` when no source tree is present, and the three
affected tests call `testthat::skip_if()` on that condition. They continue to
run, and to guard, when the suite is run from a source checkout.

This was verified by reproducing the failing condition: installing the package
into a clean library and running the suite from an isolated directory with no
source tree reachable. The three audits are reported as skipped and the suite
passes.

## Other changes in this version

**One new compiled routine.** `src/corr_gram.cpp` adds
`fgee_gram_axis_sums()`, a single pass over the cluster design forming one axis
of sums. It is registered through `RcppExports` in the usual way; the
registration audit above was extended to cover it. It allocates one
`NumericMatrix` and contains no pointer arithmetic beyond indexing into that
matrix and the input.

It supports a faster route to the working statistics when one working
correlation axis is exchangeable and the other independent, using the
Sherman-Morrison structure of the exchangeable inverse. The operator is
mathematically identical to the existing path. Agreement with the previous
implementation was checked across 96 configurations spanning positive, zero and
negative correlations up to the positive-definiteness boundary, with a worst
relative difference of 1.3e-15. Sixteen inadmissible structures were confirmed
to fall through to the previous code unchanged, and the route can be disabled
entirely with `options(fastFGEE.corr.gram = FALSE)`.

**Vignette.** The walkthrough for starting from a user-supplied
`refund::pffr()` fit again covers irregular functional grids, and a new
subsection discusses structured working correlations (stationary AR(p) and
Matern models on regular grids, and FPCA) as possible future extensions.

## Test environments

* local: x86_64 Linux, R 4.4.3 -- `R CMD check --as-cran`
* installed-package test run reproducing the CRAN flavor condition (above)
* r-universe: Windows (R 4.5.3, 4.6.1, 4.7.0), macOS (R 4.5.3, 4.6.1),
  Linux (R 4.6.1, 4.7.0)

## R CMD check results

0 errors | 0 warnings | 3 notes

* "unable to verify current time" -- the local check machine has no reachable
  time service.
* "Compilation used the following non-portable flag(s): '-march=nocona'" --
  from the local conda toolchain's default CXXFLAGS. `src/Makevars` sets only
  `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` and specifies no
  architecture flags. This note did not appear on the CRAN pretest machines.
* "Days since last update" may appear; this submission fixes the check failures
  reported for 0.2.0.

## Downstream dependencies

There are no reverse dependencies on CRAN.
