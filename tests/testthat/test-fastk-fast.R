test_that("fast gradient fastK reaches the same criterion as the stock selector", {
  skip_if_not_installed("mgcv")
  for (fam_nm in c("binomial", "poisson", "gamma")) {
    fam <- switch(fam_nm, binomial = stats::binomial(), poisson = stats::poisson(),
                  gamma = stats::Gamma(link = "log"))
    fx <- make_working_fixture(N = 8L, n_long = 5L, n_fun = 8L, p = 4L)
    fx$data[, Y := switch(fam_nm,
      binomial = stats::rbinom(.N, 1L, 0.4),
      poisson  = stats::rpois(.N, 1.2),
      gamma    = stats::rgamma(.N, shape = 3, scale = 0.4))]
    fit <- make_fake_fit(p = 4L, nsdf = 1L, q = 1L, family = fam)
    st <- fgee_build_working_stats(fx$data, fx$namesd, retain = "scores",
                                   corr_fn = "ar1", corr_long = "exchangeable",
                                   copy_dt = FALSE)
    prep <- fgee_fastk_prepare(st, fx$data, fx$namesd, fit, K = 3L, seed = 1L,
                               exact = FALSE, memory = "balanced")
    slow <- fgee_tune_fastk_grad(prep, verbose = FALSE)
    fast <- fgee_tune_fastk_grad_fast(prep, verbose = FALSE)
    # Both minimise the same criterion, so the fast path must not end up at a
    # materially worse objective. Lambda itself is only weakly identified on the
    # flat parts of the surface and is deliberately not compared.
    expect_lt(fast$score, slow$score * (1 + 1e-4) + 1e-8)
    expect_length(fast$lambda, prep$q)
    expect_true(all(is.finite(fast$lambda)) && all(fast$lambda > 0))
  }
})

test_that("the compiled kernel reproduces the R objective exactly", {
  skip_if_not(fgee_fastk_kernel_ok(), "compiled kernel unavailable")
  fx <- make_working_fixture(N = 7L, n_long = 4L, n_fun = 9L, p = 4L)
  fx$data[, Y := stats::rbinom(.N, 1L, 0.45)]
  fit <- make_fake_fit(p = 4L, nsdf = 1L, q = 1L, family = stats::binomial())
  st <- fgee_build_working_stats(fx$data, fx$namesd, retain = "scores",
                                 corr_fn = "ar1", corr_long = "ar1",
                                 copy_dt = FALSE)
  prep <- fgee_fastk_prepare(st, fx$data, fx$namesd, fit, K = 3L, seed = 1L,
                             exact = FALSE, memory = "balanced")
  set.seed(4)
  for (i in 1:5) {
    lam <- prep$lambda_frem * 10^stats::runif(prep$q, -3, 3)
    a <- fgee_fastk_score_grad(prep, lam, need_gradient = TRUE)
    b <- fastFGEE:::.fgee_fastk_score_grad_kernel(prep, lam, need_gradient = TRUE)
    expect_equal(b$score, a$score, tolerance = 1e-12)
    expect_equal(b$gradient, a$gradient, tolerance = 1e-9)
  }
})

test_that("disabling the kernel cannot change the selected smoothing parameters", {
  skip_if_not(fgee_fastk_kernel_ok(), "compiled kernel unavailable")
  fx <- make_working_fixture(N = 8L, n_long = 4L, n_fun = 8L, p = 4L)
  fx$data[, Y := stats::rpois(.N, 1.1)]
  fit <- make_fake_fit(p = 4L, nsdf = 1L, q = 1L, family = stats::poisson())
  st <- fgee_build_working_stats(fx$data, fx$namesd, retain = "scores",
                                 corr_fn = "ar1", corr_long = "exchangeable",
                                 copy_dt = FALSE)
  prep <- fgee_fastk_prepare(st, fx$data, fx$namesd, fit, K = 3L, seed = 1L,
                             exact = FALSE, memory = "balanced")
  with_k <- fgee_tune_fastk_grad_fast(prep, kernel = TRUE, verbose = FALSE)
  no_k <- fgee_tune_fastk_grad_fast(prep, kernel = FALSE, verbose = FALSE)
  expect_equal(with_k$lambda, no_k$lambda, tolerance = 1e-10)
  expect_equal(with_k$score, no_k$score, tolerance = 1e-12)
})

test_that("stage 1 reports whether its optimum is interior", {
  fx <- make_working_fixture(N = 6L, n_long = 4L, n_fun = 7L, p = 4L)
  fx$data[, Y := stats::rbinom(.N, 1L, 0.5)]
  fit <- make_fake_fit(p = 4L, nsdf = 1L, q = 1L, family = stats::binomial())
  st <- fgee_build_working_stats(fx$data, fx$namesd, retain = "scores",
                                 corr_fn = "ar1", corr_long = "ar1",
                                 copy_dt = FALSE)
  prep <- fgee_fastk_prepare(st, fx$data, fx$namesd, fit, K = 3L, seed = 1L,
                             exact = FALSE, memory = "balanced")
  res <- fgee_tune_fastk_grad_fast(prep, verbose = FALSE)
  # A non-interior stage-1 optimum means the ray search was truncated rather
  # than converged, so this flag must always be present and TRUE here.
  expect_true(isTRUE(res$stage1_interior))
  expect_true(res$stage2_mode %in% c("cartesian", "coordinate"))
})
