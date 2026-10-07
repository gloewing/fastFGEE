# fastFGEE 0.3.0.9006 development notes

This tree was reconstructed from the validated `0.3.0.9004` source rather than
repairing the failed `0.3.0.9005` tree in place. Only the independently audited
changes described below were ported forward.

## Public estimator contract

The exported `fgee()` function fits the one-step estimator only. The historical
exact Gaussian, pffr-only, legacy-engine, and fully iterated implementations are
retained in the unexported `.fgee_fit_internal()` development interface so that
legacy-versus-optimized tests and future methodological experiments remain
possible. The exported signature does not contain `exact`, `gee.fit`,
`max.iter`, `tune.method`, or `working.engine`, and attempts to pass these names
through `...` fail explicitly.

## Core working-statistics identity

For cluster `i`, let

```
r_i = (y_i - mu_i) / sqrt(v_i)
Z_i = diag(muprime_i / sqrt(v_i)) X_i.
```

The optimized engine forms

```
W_i = t(Z_i) R_i^{-1} Z_i
d_i = t(Z_i) R_i^{-1} r_i
```

by applying `R_i^{-1}` once to `cbind(Z_i, r_i)`. The public fit always uses
one coefficient update and the optimized engine.

## Nuisance integration

The typed Gamma-dispersion, beta-precision, and NB2-size implementation from
`0.3.0.9004` is preserved. The former `zzz_nuisance_integration.R` load-order
capture has been removed. Core implementations now have explicit internal
names and the active wrappers call those names directly. Nuisance parameters
remain fixed during smoothing-parameter tuning and are updated only at the
designated final working-state calculation.

## Numerical code

`solve_pd()` now uses an internal Rcpp/LAPACK implementation of Cholesky solve
and inverse operations (`dpotrf`, `dpotrs`, and `dpotri`) with a base-R fallback.
The package no longer calls `sanic`.

Irregularly sampled continuous-time AR(1) working correlations use the exact
tridiagonal Markov precision for

```
Corr{Z(t_j), Z(t_k)} = rho^abs(t_j - t_k),  0 <= rho < 1.
```

The precision can be applied to multiple right-hand sides in linear time. The
implementation is independent and follows the formulas in Allevius (2018); no
source from the archived `irregulAR1` package was copied. The package no longer
calls or declares `irregulAR1`. A bounded profile-likelihood estimator is kept
as an internal method-development helper; the validated public correlation
update schedule is otherwise unchanged.

All Rcpp entry points are listed in both generated registration files. Rcpp is
a direct `Imports` and `LinkingTo` dependency. The stale `sourceCpp` import was
replaced by the standard `evalCpp` namespace import.

## Smoothing methods and memory layouts

The automatic dispatch is unchanged:

* Gaussian identity: `fastk_staged`.
* Every supported non-Gaussian family, including NB2 and beta: `fastk_grad_fast`.

`qreml_fastk` is retained for reproducibility but normally reaches the same
fastK solution by a slower route. `sandwich_qreml` remains an opt-in
experimental selector because its performance can deteriorate under working
correlation misspecification.

`fastk.memory = "lowmem"` is deprecated. Fresh-process measurements found no
regime in which it reduced peak RSS, and it cannot use the compiled fastK
kernel because only the `balanced` layout constructs `prep$X_eval`.

## Inference wording

Documentation now distinguishes pointwise intervals from optional simultaneous
bands. The wild maximum-statistic band is a simultaneous-inference procedure,
not a promise of nominal finite-sample coverage. Reliability may be poor with
few or high-leverage independent clusters or unstable nuisance/correlation
estimates.

## Internal initial-fit baseline scale

The internal initial-fit path (`gee.fit = FALSE`) holds the supplied `pffr()`
coefficients fixed and builds sandwich/bootstrap quantities from fastFGEE's
working scores. mgcv's `sp` penalise the total deviance/RSS; the one-step bread
`Wbar + P` is on the averaged cluster scale, so the penalty is
`S_lambda / (N * c_work)` with `c_work` the dispersion actually multiplied into
the working variance.

`c_work` is not always `fit$sig2`. It is `fit$sig2` for Gaussian (grid-specific
working variances, retained); the working dispersion for Gamma; and `1` for
binomial, Poisson, the quasi families on the current fixed-nuisance path,
negative binomial and beta. The quasi case is the discriminating one: its
working variance omits any extra scalar, so the divisor is `N`, not
`N * fit$sig2`.

Two engine-level facts worth remembering when writing tests here:

* `working0` is populated by the optimized engine; the legacy engine returns the
  per-cluster lists `wi0` / `di0` instead, so a test that needs the bread must
  read whichever the engine actually filled.
* For Gamma the two engines may use different working dispersions. The legacy
  path estimates one from the current fitted mean (and warns), rather than
  taking `fit$sig2`. Each engine's returned penalty is consistent with its own
  dispersion, so do not assert equality of the two engines' covariances.

Regression coverage: `tests/testthat/test-initial-penalty-scale-helper.R`
(dispatch, runs in ordinary checks) and
`tests/testthat/test-initial-penalty-scale-families.R` (Gamma and quasi-Poisson
end-to-end against an independently assembled penalty and centred cluster
sandwich). Scope limits are recorded in `man/fastFGEE-initial-fit.Rd`: this path
is not universally the exact sandwich of the original pffr likelihood estimator,
and the Gaussian exact-weighting and beta-likelihood-score questions are
tracked separately.

## Validation policy

The prior `0.3.0.9004` numerical and family validation remains the baseline.
New acceptance work for this port consists of:

1. Rcpp registration and package-build checks;
2. dense-reference tests for SPD solves and irregular AR(1) precision;
3. source audits proving that `sanic` and `irregulAR1` executable calls are gone;
4. public-signature tests proving that unsupported estimator controls are not
   exported;
5. the existing full test suite and `R CMD check` in an R-enabled environment.

The current construction environment did not provide R or Rscript. Therefore,
source-level and independent numerical validation can be completed here, but a
real `R CMD build`, installation, testthat run, and `R CMD check` must be run on
the supplied R validation script before this candidate is promoted.
