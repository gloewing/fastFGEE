test_that("optimized one-step engine reproduces the legacy engine", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("SuperGauss")

  data("d", package = "fastFGEE")
  fit0 <- refund::pffr(
    Y ~ X1 + X2,
    data = d,
    family = binomial(link = "logit"),
    algorithm = "bam",
    discrete = TRUE,
    bs.yindex = list(bs = "ps", k = 7, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 7, m = c(2, 1))
  )
  grids <- list(
    c(0.1, 1, 10),
    c(0.1, 1, 10),
    c(0.5, 1, 2)
  )

  legacy <- fastFGEE:::.fgee_fit_internal(
    Y ~ X1 + X2,
    data = d,
    cluster = "ID",
    family = binomial(link = "logit"),
    time = "time",
    corr_long = "exchangeable",
    corr_fn = "independent",
    pffr.mod = fit0,
    cv.grid = grids,
    joint.CI = FALSE,
    working.engine = "legacy",
    sp.method = "legacy",
    rho.smooth = FALSE
  )

  optimized <- fastFGEE::fgee(
    Y ~ X1 + X2,
    data = d,
    cluster = "ID",
    family = binomial(link = "logit"),
    time = "time",
    corr_long = "exchangeable",
    corr_fn = "independent",
    pffr.mod = fit0,
    cv.grid = grids,
    joint.CI = FALSE,
    sp.method = "fastk_staged",
    working.retain = "full",
    corr.solver = "supergauss",
    rho.smooth = FALSE,
    verbose.tuning = FALSE
  )

  expect_equal(optimized$lambda, as.numeric(legacy$lambda), tolerance = 1e-8)
  expect_equal(optimized$beta, legacy$beta, tolerance = 2e-7)
  expect_equal(optimized$vb, legacy$vb, tolerance = 5e-7)
  # The fit above is made as `fastFGEE::fgee(...)`, so match.call() records
  # `::`(fastFGEE, fgee) and as.character() yields three elements. Assert the
  # function name, which is what this line exists to check: a regression
  # would record an anonymous closure instead.
  expect_identical(tail(as.character(optimized$call[[1L]]), 1L), "fgee")
  expect_equal(optimized$n_iter, 1L)
  expect_false(isTRUE(optimized$exact))
  expect_true(is.list(optimized$nuisance))
})

test_that("optimized model update restores the original pffr row order", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")

  data("d", package = "fastFGEE")
  set.seed(901)
  ds <- d[sample(seq_len(nrow(d))), , drop = FALSE]
  fit0 <- refund::pffr(
    Y ~ X1 + X2,
    data = ds,
    family = binomial(link = "logit"),
    algorithm = "bam",
    discrete = TRUE,
    bs.yindex = list(bs = "ps", k = 7, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 7, m = c(2, 1))
  )
  grids <- list(c(0.1, 1, 10), c(0.1, 1, 10), c(0.5, 1, 2))
  ans <- fastFGEE::fgee(
    Y ~ X1 + X2,
    data = ds,
    cluster = "ID",
    family = binomial(link = "logit"),
    time = "time",
    corr_long = "exchangeable",
    corr_fn = "independent",
    pffr.mod = fit0,
    cv.grid = grids,
    joint.CI = FALSE,
    sp.method = "fastk_staged",
    working.retain = "scores",
    rho.smooth = FALSE,
    verbose.tuning = FALSE
  )

  X <- suppressWarnings(stats::model.matrix(fit0))
  eta <- as.numeric(X %*% ans$beta)
  expect_equal(ans$model$linear.predictors, eta, tolerance = 1e-10)
  expect_equal(ans$model$fitted.values, stats::plogis(eta), tolerance = 1e-10)
})
