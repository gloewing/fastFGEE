# .fgee_initial_penalty_scale(): family dispatch for the gee.fit = FALSE path ---
#
# mgcv's smoothing parameters penalise the TOTAL deviance/RSS, while the one-step
# bread is on the averaged cluster scale Wbar + P.  The helper returns the factor
# G * c_work that converts between them, where c_work is the dispersion actually
# multiplied into fastFGEE's working variance -- not necessarily fit$sig2.
#
# These are pure dispatch tests on the helper contract (family + sig2 + nuisance),
# so they are cheap and run in ordinary checks.  The fitted-model behaviour that
# establishes which nuisance record and working variance the engines really use is
# tested in test-initial-penalty-scale-families.R.

fake_fit <- function(family, sig2 = 1) list(family = family, sig2 = sig2)

test_that("Gamma scale is G times the dispersion supplied by the working nuisance", {
  for (phi in c(0.4, 2.5)) {
    for (G in c(17, 40)) {
      expect_equal(
        fastFGEE:::.fgee_initial_penalty_scale(
          fake_fit(stats::Gamma(link = "log"), sig2 = phi), G,
          list(family = "gamma", dispersion = phi)
        ),
        G * phi,
        label = paste("Gamma phi =", phi, "G =", G)
      )
    }
  }
})

test_that("the Gamma scale follows the working nuisance, not fit$sig2", {
  # Dispatch guard: the documented source of c_work is the dispersion actually
  # multiplied into the working variance, carried on the nuisance record.  If an
  # implementation silently read fit$sig2 instead, this case would return 30 * 1.9.
  #
  # This is a deliberately constructed contract case.  In the real fitted Gamma
  # fixture the two agree (see test-initial-penalty-scale-families.R); nothing
  # here asserts an observed disagreement in practice.
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(
      fake_fit(stats::Gamma(link = "log"), sig2 = 1.9), 30,
      list(family = "gamma", dispersion = 0.4)
    ),
    30 * 0.4
  )
})

test_that("the quasi-Poisson fixed-nuisance path scales by G alone", {
  # The current quasi path does NOT multiply an extra scalar dispersion into the
  # working variance, so c_work = 1 and the factor is G even when the fitted sig2
  # is far from one.  A mistaken G * fit$sig2 implementation would return 25 * 1.7.
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(
      fake_fit(stats::quasipoisson(link = "log"), sig2 = 1.7), 25, NULL
    ),
    25
  )
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(
      fake_fit(stats::quasibinomial(link = "logit"), sig2 = 2.3), 12, NULL
    ),
    12
  )
})

test_that("Gaussian keeps the G * fit$sig2 convention", {
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(
      fake_fit(stats::gaussian(), sig2 = 0.64), 30, NULL
    ),
    30 * 0.64
  )
})

test_that("negbinomial and beta take no extra size/precision factor", {
  # theta and the beta precision already enter V(mu) exactly as in mgcv's
  # extended families, so no further global division is applied.
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(
      fake_fit(mgcv::nb(link = "log"), sig2 = 1), 22,
      list(family = "negbinomial", value = 3)
    ),
    22
  )
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(
      fake_fit(mgcv::betar(link = "logit"), sig2 = 1), 22,
      list(family = "beta", value = 7, dispersion = 2.5)
    ),
    22
  )
})

test_that("canonical one-parameter families scale by G", {
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(fake_fit(stats::poisson()), 31, NULL), 31
  )
  expect_equal(
    fastFGEE:::.fgee_initial_penalty_scale(fake_fit(stats::binomial()), 31, NULL), 31
  )
})

test_that("N must be a positive cluster count", {
  expect_error(
    fastFGEE:::.fgee_initial_penalty_scale(fake_fit(stats::poisson()), 0, NULL),
    "positive number of clusters"
  )
})
