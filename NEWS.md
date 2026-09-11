# fastFGEE 0.2.0

First release since 0.1.0. The public interface is largely unchanged --
`fgee()` and `fgee.plot()` keep their meaning and their required arguments --
but the estimation engine, the working-correlation machinery and the dependency
footprint have all changed substantially.

## Breaking changes

* Several arguments that selected non-default estimation workflows have been
  removed from `fgee()`: `exact`, `gee.fit`, `max.iter` and `tune.method`.
  Passing any of them (including through `...`) is now an error rather than
  being silently ignored, so existing scripts fail loudly instead of quietly
  changing meaning. `fgee()` fits the one-step estimator, which is the method
  the package documents and the accompanying paper analyses.

* `'irregulAR1'` is no longer used. 0.1.0 could optionally call this archived
  CRAN package for irregularly spaced AR(1) precision matrices; the operation
  is now implemented inside the package, so no archived dependency is needed.
  Note that the internal routine returns the precision matrix directly, whereas
  `irregulAR1` returned it scaled by `(1 - rho^2)`.

* `'sanic'` is no longer used; symmetric positive-definite solves now go
  through LAPACK directly.

## New features

* New `rho.pool` argument to `fgee()` controls whether the surviving
  working-correlation parameter is pooled to a cluster-level scalar when only
  one direction is correlated. The default `"fn"` estimates a single pooled
  functional correlation when `corr_long = "independent"`, instead of one value
  per longitudinal observation. This restores the separable (Kronecker)
  structure, so the correlation inverse is applied once per cluster rather than
  once per longitudinal observation. `rho.pool = "none"` reproduces the previous
  behaviour. The longitudinal direction is deliberately not pooled by default.

* Compiled kernels for the performance-critical operations: symmetric
  positive-definite Cholesky solves, exact tridiagonal precision operations for
  irregularly sampled continuous-time AR(1) working correlations, and a batched
  Kronecker inverse applied to the working statistics. The package now has a
  `src/` directory and requires compilation.

* Smoothing-parameter selection gained several selectors, exposed through
  `sp.method`: staged and analytic-gradient fast cluster cross-validation, and
  an experimental sandwich-scaled working restricted quasi-likelihood selector.
  The default `"auto"` dispatches on the family.

* New arguments give control over tuning and memory behaviour: `fastk.K`,
  `fastk.seed`, `fastk.memory`, `fastk.start`, `fastk.kernel`, `working.retain`,
  `corr.solver`, `qreml.phi.method`, `qreml.phi.fixed`, `qreml.phi.weight`,
  `qreml.phi.clip`, `keep.tuning.workspace`, `keep.data`, `keep.initial.fit`,
  `keep.working.stats` and `verbose.tuning`.

* Working correlations in the functional direction now include `"fpca"`.

* New `print()` method for `fgee_working_stats` objects.

## Improvements

* Working statistics for each cluster are accumulated in a single pass, which
  removes the repeated construction of large intermediate matrices that
  dominated cost in 0.1.0 for large cluster sizes.

* Negative binomial, Gamma and beta nuisance parameters are estimated
  consistently across the tuning and inference paths.

* Fitted objects no longer inline the model frame into `fit$call`, so saved
  fits are substantially smaller.

* Removed explicit `gc()` calls from the estimation code. These forced a full
  garbage collection after each working-statistics build, costing roughly 18%
  of fit time on the legacy engine with no memory benefit; results are
  unchanged.

## Dependency changes

* `Rcpp` moves from `Suggests` to `Imports`, and is added to `LinkingTo`.
  `NeedsCompilation` is now `yes`: a C++ toolchain is required to install from
  source. There is no pure-R install-time fallback.
* `Rfast` is no longer a dependency.
* `SuperGauss` moves from `Imports` to `Suggests`; it is used only when
  `corr.solver = "supergauss"` is requested.
* `irregulAR1` and `sanic` are removed from `Suggests`.
* `testthat` (>= 3.0.0) added to `Suggests`; the package now ships a test suite.

## Documentation

* Added a note to `corr_fn`/`corr_long` that setting exactly one direction to
  `"independent"` forfeits the separable structure, so a simpler working
  correlation is not necessarily cheaper to fit.
* `NEWS.md` added.
