test_that("qREML analytic gradient matches finite differences", {
  fx <- make_working_fixture(N = 30L, n_long = 2L, n_fun = 4L, p = 7L, seed = 91)
  fit <- make_fake_fit(p = 7L, q = 2L, family = gaussian())
  st <- fastFGEE:::fgee_build_working_stats(fx$data, fx$namesd, retain = "scores")
  prep <- fastFGEE:::fgee_qreml_prepare(st, fit)
  info <- fastFGEE:::fgee_qreml_information(prep)
  phi <- min(max(info$phi_pen, 1), 8)
  ne <- st$N / phi
  u <- log10(c(0.9, 1.4) / prep$lambda_frem)
  analytic <- fastFGEE:::fgee_qreml_score_grad(
    prep, prep$lambda_frem * 10^u, ne, TRUE
  )$gradient
  numeric <- finite_diff(function(v) {
    fastFGEE:::fgee_qreml_score_grad(
      prep, prep$lambda_frem * 10^v, ne, FALSE
    )$score
  }, u)
  expect_equal(analytic, numeric, tolerance = 2e-5)
})

test_that("one-sided sandwich scaling never claims extra information", {
  fx <- make_working_fixture(N = 20L, p = 7L, seed = 100)
  fit <- make_fake_fit(p = 7L, q = 2L)
  st <- fastFGEE:::fgee_build_working_stats(fx$data, fx$namesd, retain = "scores")
  prep <- fastFGEE:::fgee_qreml_prepare(st, fit)
  ans <- fastFGEE:::fgee_tune_qreml(
    prep, phi_method = "fixed", phi_fixed = 0.4,
    phi_clip = c(1, 8), coarse = c(0), n_starts = 1L, maxit = 2L
  )
  expect_equal(ans$phi_used, 1)
  expect_equal(ans$n_eff, st$N)
})

test_that("exact Gaussian qREML uses the exact-estimator penalty scale", {
  fx <- make_working_fixture(N = 24L, n_long = 3L, n_fun = 4L, p = 7L,
                             seed = 191)
  fit <- make_fake_fit(p = 7L, q = 2L, family = gaussian())
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores", exact_gaussian = TRUE,
    beta0 = fit$coefficients
  )
  prep <- fastFGEE:::fgee_qreml_prepare(st, fit, exact = TRUE)
  lambda <- c(0.7, 1.8)
  ev <- fastFGEE:::fgee_qreml_score_grad(
    prep, lambda, n_eff = st$N, need_gradient = FALSE
  )
  P <- fastFGEE:::penalty_from_setup(prep$ps, lambda)
  exact_beta <- fastFGEE:::solve_pd(st$W_sum + P, st$d_exact_sum)

  expect_equal(prep$penalty_scale, 1 / st$N)
  expect_equal(ev$beta, as.numeric(exact_beta), tolerance = 1e-11)

  tuned <- fastFGEE:::fgee_tune_qreml(
    prep, phi_method = "fixed", phi_fixed = 1,
    coarse = c(0), n_starts = 1L, maxit = 2L
  )
  expect_equal(
    tuned$penalty_mat,
    fastFGEE:::penalty_from_setup(prep$ps, tuned$lambda),
    tolerance = 1e-12
  )
  expect_equal(
    tuned$penalty_effective,
    tuned$penalty_mat / st$N,
    tolerance = 1e-12
  )
})

test_that("qREML uses aggregate statistics and does not require retained scores", {
  fx <- make_working_fixture(N = 30L, n_long = 2L, n_fun = 4L, p = 7L,
                             seed = 216)
  fit <- make_fake_fit(p = 7L, q = 2L, family = gaussian())
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "aggregate"
  )
  expect_null(st$D)
  prep <- fastFGEE:::fgee_qreml_prepare(st, fit)
  info <- fastFGEE:::fgee_qreml_information(prep)
  expect_true(is.finite(info$phi_pen))
  ans <- fastFGEE:::fgee_tune_qreml(
    prep, coarse = c(-1, 0, 1), n_starts = 1L, maxit = 3L,
    verbose = FALSE
  )
  expect_true(all(is.finite(ans$lambda)))
  expect_true(all(ans$lambda > 0))
})
