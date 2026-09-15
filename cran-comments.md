## Resubmission

This is a resubmission. The previous submission was rejected by the incoming
pretest with two NOTEs; both are addressed below.

### 1. "no visible global function definition for 'head', 'setNames', 'tail'"

Fixed. `utils::head()`, `utils::tail()` and `stats::setNames()` were used
without being imported. The following have been added to NAMESPACE:

    importFrom(stats, setNames)
    importFrom(utils, head)
    importFrom(utils, tail)

`R CMD check` now reports no undefined global functions or variables.

This was not caught before submission because the local check machine was
missing the recommended package 'codetools', which caused
"checking R code for possible problems" to be skipped silently. 'codetools'
has been installed and the check re-run.

### 2. "Possibly misspelled words in DESCRIPTION: tridiagonal"

This is a false positive. "Tridiagonal" is standard linear-algebra
terminology for a matrix whose non-zero entries lie on the main diagonal and
the two adjacent diagonals. It is used correctly in the Description to
describe the exact precision-matrix operations for continuous-time AR(1)
working correlations. No change has been made.

## Test environments

* local: x86_64 Linux, R 4.4.3 -- `R CMD check --as-cran`
* win-builder incoming pretest: Windows and Debian (previous submission) --
  both installed, loaded, and passed tests, examples, vignette rebuild and
  PDF manual
* r-universe: Windows (R 4.5.3, 4.6.1, 4.7.0), macOS (R 4.5.3, 4.6.1),
  Linux (R 4.6.1, 4.7.0), wasm (R 4.6.0) -- all build and install cleanly

## R CMD check results

0 errors | 0 warnings | 2 notes locally

Both remaining local notes are properties of the check machine, not the
package, and did not appear on the CRAN pretest machines:

* "unable to verify current time" -- the machine has no reachable time service.
* "Compilation used the following non-portable flag(s): '-march=nocona'" --
  this comes from the local conda toolchain's default CXXFLAGS.
  `src/Makevars` sets only
  `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` and specifies no
  architecture flags.

## Changes in this version (0.2.0), first release since 0.1.0

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
