# The frozen scalar nuisance on the fastK workspace.
#
# Two properties are asserted here. First, that the criterion and the working
# statistics share one number -- they must, because Bhat[, k] is built from
# prep$Wbar and prep$D, which are functions of the working variance, so scoring
# those coefficients under a different nuisance would not be a held-out
# predictive loss for the model being fitted. Second, that tuning never moves
# it: that is the assertion which would catch a future refactor re-estimating
# the nuisance per candidate smoothing parameter.

gamma_prep <- function(nuisance_value = 0.4, ..., seed = 481) {
  fx <- make_working_fixture(N = 8L, n_long = 4L, n_fun = 6L, p = 5L,
                             seed = seed)
  rec <- structure(
    list(family = "gamma", parameter = "dispersion",
         value = nuisance_value, dispersion = nuisance_value,
         method = "fixed", source = "test"),
    class = c("fgee_nuisance", "list")
  )
  data.table::setattr(fx$data, "nuisance", rec)
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "scores", copy_dt = TRUE
  )
  fit <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L,
                       family = stats::Gamma(link = "log"))
  list(prep = fastFGEE:::fgee_fastk_prepare(st, fx$data, fx$namesd, fit,
                                            K = 4L, seed = 1L, ...),
       working = st, fit = fit, fx = fx)
}

test_that("the workspace inherits the nuisance that built the working statistics", {
  g <- gamma_prep(nuisance_value = 0.4)
  expect_identical(g$prep$nuisance_source, "working$nuisance")
  expect_equal(g$prep$nuisance_value, 0.4)
  # the number the criterion holds is the number the working statistics used
  expect_equal(g$prep$nuisance_value,
               fastFGEE:::.fgee_nuisance_value(g$working, "dispersion", NA_real_))
})

test_that("an explicit override is recorded as such", {
  g <- gamma_prep(nuisance_value = 0.4, nuisance = 7.5)
  expect_identical(g$prep$nuisance_source, "override")
  expect_equal(g$prep$nuisance_value, 7.5)
})

test_that("families without a scalar nuisance carry none", {
  for (fam in list(stats::gaussian(), stats::binomial(), stats::poisson())) {
    fx <- make_working_fixture(N = 8L, n_long = 4L, n_fun = 6L, p = 5L)
    st <- fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
      retain = "scores", copy_dt = TRUE)
    fit <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L, family = fam)
    prep <- fastFGEE:::fgee_fastk_prepare(st, fx$data, fx$namesd, fit,
                                          K = 4L, seed = 1L)
    expect_true(is.na(prep$nuisance_value), info = fam$family)
    expect_true(is.na(prep$nuisance_source), info = fam$family)
  }
})

test_that("the gamma criterion is inert to the nuisance, so the channel is inert", {
  # Gamma dispersion cancels out of the fastK criterion: it enters the
  # log-likelihood multiplicatively plus eta-free terms. So plumbing it through
  # must not move the gamma score by even one bit -- this is the regression
  # guard for step 2 of the change.
  lam <- c(0.7, 3.2)
  a <- gamma_prep(nuisance_value = 0.4, nuisance = 1e-3)
  b <- gamma_prep(nuisance_value = 0.4, nuisance = 1e3)
  sa <- fastFGEE:::fgee_fastk_score_grad(a$prep, lam, need_gradient = TRUE)
  sb <- fastFGEE:::fgee_fastk_score_grad(b$prep, lam, need_gradient = TRUE)
  expect_identical(sa$score, sb$score)
  expect_identical(sa$gradient, sb$gradient)
})

test_that("tuning does not move the frozen nuisance", {
  g <- gamma_prep(nuisance_value = 0.4)
  before_value <- g$prep$nuisance_value
  before_record <- g$prep$working$nuisance
  invisible(fastFGEE:::fgee_tune_fastk_grad_fast(g$prep, verbose = FALSE))
  expect_identical(g$prep$nuisance_value, before_value)
  expect_identical(g$prep$working$nuisance, before_record)
})
