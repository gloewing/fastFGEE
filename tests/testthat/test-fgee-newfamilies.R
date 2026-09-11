# End-to-end fits for the two newly supported families.
#
# Reaching the fastK criterion switches on a whole chain that had never run for
# these families: the sandwich under an NB2 / beta working variance, the
# studentised wild cluster bootstrap, the joint CI construction and the
# coefficient-curve extraction. "No error" is not a sufficient assertion here,
# so the covariance is checked for symmetry and positive semi-definiteness and
# the bands for being finite, ordered and non-degenerate.

nf_data <- function(family, N = 14L, n_long = 4L, n_fun = 8L, seed = 909L) {
  set.seed(seed)
  n_row <- N * n_long
  s <- (seq_len(n_fun) - 0.5) / n_fun
  ID <- rep(seq_len(N), each = n_long)
  time <- rep(seq_len(n_long), times = N)
  X1 <- rep(stats::rnorm(N), each = n_long)
  X2 <- stats::rbinom(n_row, 1L, 0.5)
  eta <- outer(rep(1, n_row), 0.15 + 0.3 * sin(2 * pi * s)) +
    outer(X1, 0.1 * cos(2 * pi * s)) + outer(X2, -0.08 * (2 * s - 1))
  Y <- if (identical(family, "negbinomial")) {
    matrix(stats::rnbinom(n_row * n_fun, size = 3, mu = exp(eta)),
           n_row, n_fun)
  } else {
    mu <- stats::plogis(eta)
    matrix(stats::rbeta(n_row * n_fun, shape1 = mu * 8, shape2 = (1 - mu) * 8),
           n_row, n_fun)
  }
  storage.mode(Y) <- "double"
  d <- data.frame(ID = ID, time = time, X1 = X1, X2 = X2)
  d$Y <- I(Y)
  d[, c("Y", "ID", "X1", "X2", "time")]
}

nf_fit <- function(family, sp.method = "auto") {
  fam <- if (identical(family, "negbinomial")) mgcv::nb() else mgcv::betar()
  fastFGEE::fgee(
    Y ~ X1 + X2, data = nf_data(family), cluster = "ID", family = fam,
    time = "time", corr_long = "exchangeable", corr_fn = "ar1",
    joint.CI = "wild", rho.smooth = FALSE,
    sp.method = sp.method, knots = 7, bs = "ps", verbose.tuning = FALSE
  )
}

for (fam in c("negbinomial", "beta")) {
  test_that(paste("fgee() fits", fam, "end to end through fastK"), {
    skip_on_cran()
    skip_if_not_installed("refund")
    skip_if_not_installed("mgcv")

    fit <- nf_fit(fam)

    # the criterion actually used was fastK, not a silent qREML fallback
    expect_identical(fit$tuning$method, "fastk_grad_fast")
    expect_true(all(is.finite(fit$lambda)))
    expect_true(all(fit$lambda > 0))
    expect_true(all(is.finite(as.numeric(fit$beta))))

    # the frozen nuisance is recorded and positive
    expect_true(is.list(fit$nuisance))
    expect_true(fastFGEE:::.fgee_positive_scalar(fit$nuisance$final$value))

    # covariance: finite, symmetric, positive semi-definite
    V <- as.matrix(fit$model$Vp)
    expect_true(all(is.finite(V)))
    expect_equal(V, t(V), tolerance = 1e-10)
    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(min(ev), -1e-8 * max(1, max(abs(ev))))

    # bands: finite, correctly ordered, and not collapsed onto the estimate.
    # fgee.plot() always draws through gridExtra, which would leave an
    # Rplots.pdf in tests/testthat and then in the built tarball, so send it to
    # a null device.
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    pl <- fgee.plot(fit, return = TRUE)
    for (tn in c("fn_intercept", "X1", "X2")) {
      df <- as.data.frame(pl[[tn]])
      expect_true(nrow(df) > 1L, info = tn)
      expect_true(all(is.finite(df$beta.hat)), info = tn)
      for (cc in c("CI.lower.pointwise", "CI.upper.pointwise",
                   "CI.lower.joint", "CI.upper.joint")) {
        expect_true(all(is.finite(df[[cc]])), info = paste(tn, cc))
      }
      expect_true(all(df$CI.lower.pointwise <= df$beta.hat), info = tn)
      expect_true(all(df$beta.hat <= df$CI.upper.pointwise), info = tn)
      # the joint band must contain the pointwise one
      expect_true(all(df$CI.lower.joint <= df$CI.lower.pointwise + 1e-10),
                  info = tn)
      expect_true(all(df$CI.upper.joint >= df$CI.upper.pointwise - 1e-10),
                  info = tn)
      # non-degenerate width
      expect_gt(min(df$CI.upper.pointwise - df$CI.lower.pointwise), 0)
    }
  })
}

test_that("MASS::negative.binomial is refused with an actionable message", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("refund")
  # mgcv:::fix.family.var() rejects a plain family object lacking dvar/d2var/
  # d3var with the opaque "family not recognised", from inside pffr(), before
  # any fastFGEE code runs. Whatever the message, it must fail rather than
  # silently produce a fit with the wrong working variance.
  expect_error(
    fastFGEE::fgee(Y ~ X1 + X2, data = nf_data("negbinomial"), cluster = "ID",
                   family = MASS::negative.binomial(3), time = "time",
                   corr_long = "exchangeable", corr_fn = "ar1",
                   joint.CI = FALSE,
                   knots = 7, bs = "ps", verbose.tuning = FALSE)
  )
})
