## Submission

This is a major update to fastFGEE, the first since 0.1.0. It replaces the
estimation engine, adds compiled routines, and removes two optional
dependencies. The exported interface (`fgee()`, `fgee.plot()`) is unchanged in
meaning and in its required arguments.

## Test environments

* local: x86_64 Linux, R 4.4.3 -- `R CMD check --as-cran`
* r-universe: Windows (R 4.5.3, 4.6.1, 4.7.0), macOS (R 4.5.3, 4.6.1),
  Linux (R 4.6.1, 4.7.0), wasm (R 4.6.0) -- all build and install cleanly

## R CMD check results

0 errors | 0 warnings | 2 notes

Both notes are properties of the local check machine, not the package:

* "unable to verify current time" -- the machine has no reachable time service.
* "Compilation used the following non-portable flag(s): '-march=nocona'" --
  this flag comes from the local conda toolchain's default CXXFLAGS, not from
  the package. `src/Makevars` sets only
  `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` and specifies no
  architecture flags.

## Changes that reviewers may wish to note

* The package now contains compiled code (`NeedsCompilation: yes`). Rcpp moves
  from `Suggests` to `Imports` and is added to `LinkingTo`. All compiled entry
  points are registered and `R_useDynamicSymbols` is disabled.

* The archived package 'irregulAR1' is no longer used. Exact tridiagonal
  precision operations for irregularly sampled continuous-time AR(1) working
  correlations are now implemented inside fastFGEE from the published formulas
  of Allevius (2018). No code from the archived package is used. The previous
  DESCRIPTION text directing users to the CRAN Archive has been removed.

* The 'sanic' dependency is no longer used. Symmetric positive-definite solves
  use LAPACK through registered compiled routines, with a base-R Cholesky
  fallback.

* 'Rfast' is no longer a dependency. 'SuperGauss' moves from `Imports` to
  `Suggests`; both it and 'RcppArmadillo' are reached only through
  `requireNamespace()` guards.

* Four arguments that selected non-default estimation workflows have been
  removed from `fgee()`: `exact`, `gee.fit`, `max.iter` and `tune.method`.
  They now raise an informative error rather than being silently ignored, so
  code written against 0.1.0 fails loudly instead of quietly changing meaning.
  This is documented in NEWS.md under "Breaking changes".

* A test suite (testthat, edition 3) is included for the first time: 103 test
  files, 2112 passing assertions, 0 failures.

## Downstream dependencies

There are no reverse dependencies on CRAN.
