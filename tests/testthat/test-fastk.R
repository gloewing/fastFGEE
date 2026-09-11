test_that("Gaussian sufficient-statistic score equals direct cluster loss", {
  fx <- make_working_fixture(N = 10L, n_long = 3L, n_fun = 5L, p = 7L, seed = 70)
  fit <- make_fake_fit(p = 7L, q = 2L, family = gaussian())
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores", exact_gaussian = TRUE,
    beta0 = fit$coefficients
  )
  prep <- fastFGEE:::fgee_fastk_prepare(
    st, fx$data, fx$namesd, fit, K = 5L, seed = 2L
  )
  ev <- fastFGEE:::fgee_fastk_score_grad(prep, c(0.7, 2.1), FALSE)
  B <- ev$beta_folds
  direct <- 0
  cid <- as.character(fx$data$cname_)
  for (k in seq_len(prep$K)) {
    ids <- prep$working$cluster_id[prep$folds$sets[[k]]]
    for (id in ids) {
      ii <- which(cid == id)
      e <- fx$data$Y[ii] - as.numeric(as.matrix(fx$data[ii, fx$namesd, with = FALSE]) %*% B[, k])
      direct <- direct + mean(e^2) / prep$ncl
    }
  }
  expect_equal(ev$score, direct, tolerance = 1e-11)
})

test_that("analytic fastK gradient matches finite differences", {
  fx <- make_working_fixture(N = 9L, n_long = 3L, n_fun = 4L, p = 7L, seed = 81)
  fit <- make_fake_fit(p = 7L, q = 2L, family = binomial())
  fx$data[, Y := rbinom(.N, 1L, 0.45)]
  st <- fastFGEE:::fgee_build_working_stats(fx$data, fx$namesd, retain = "scores")
  prep <- fastFGEE:::fgee_fastk_prepare(
    st, fx$data, fx$namesd, fit, K = 3L, seed = 4L, memory = "balanced"
  )
  u <- log10(c(0.8, 1.7) / prep$lambda_frem)
  analytic <- fastFGEE:::fgee_fastk_score_grad(
    prep, prep$lambda_frem * 10^u, TRUE
  )$gradient
  numeric <- finite_diff(function(v) {
    fastFGEE:::fgee_fastk_score_grad(
      prep, prep$lambda_frem * 10^v, FALSE
    )$score
  }, u)
  expect_equal(analytic, numeric, tolerance = 2e-5)
})

test_that("Gamma family names normalize correctly", {
  expect_equal(fastFGEE:::.fgee_family_key(Gamma(link = "log")), "gamma")
  expect_equal(fastFGEE:::.fgee_family_key("Gamma"), "gamma")
})

test_that("compact staged fastK reproduces the legacy objective and selection", {
  for (fam in list(gaussian(), binomial())) {
    fx <- make_working_fixture(N = 12L, n_long = 3L, n_fun = 5L, p = 7L,
                               seed = if (fastFGEE:::.fgee_family_key(fam) == "gaussian") 141 else 142)
    fit <- make_fake_fit(p = 7L, q = 2L, family = fam)
    if (fastFGEE:::.fgee_family_key(fam) == "binomial") {
      fx$data[, Y := rbinom(.N, 1L, 0.45)]
    }
    st <- fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, retain = "full", corr_solver = "exact"
    )
    folds <- fastFGEE:::.fgee_make_cluster_folds(st$cluster_id, K = 4L, seed = 11L)
    cv_grid <- list(c(0.1, 1, 10), c(0.1, 1, 10), c(0.5, 1, 2))
    prep <- fastFGEE:::fgee_fastk_prepare(
      st, fx$data, fx$namesd, fit, K = 4L, seed = 11L,
      sets = folds$sets, memory = "speed"
    )
    compact <- fastFGEE:::fgee_tune_fastk_staged(prep, cv_grid, verbose = FALSE)

    legacy_grid <- fastFGEE:::.fgee_prepare_staged_grids(cv_grid, prep$lambda_frem)
    legacy <- fastFGEE:::fun.gee1step.cv(
      w = st$W,
      d = fastFGEE:::.fgee_working_dlist(st),
      grid = legacy_grid,
      data = fx$data,
      namesd = fx$namesd,
      cname_ = "cname_",
      fit.initial = fit,
      cv = "fastkfold",
      K = 4L,
      sets = folds$sets,
      folds.list = replicate(4L, fit$coefficients),
      seed = 11L,
      exact = FALSE
    )

    expect_equal(compact$mse, legacy$mse, tolerance = 5e-10)
    expect_equal(compact$se, legacy$se, tolerance = 5e-10)
    expect_equal(unname(compact$lambda.star), unname(legacy$lambda.star), tolerance = 1e-12)
  }
})

test_that("exact Gaussian fastK uses exact cluster right-hand-side totals", {
  fx <- make_working_fixture(N = 10L, n_long = 3L, n_fun = 5L, p = 7L,
                             seed = 181)
  fit <- make_fake_fit(p = 7L, q = 2L, family = gaussian())
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores", exact_gaussian = TRUE,
    beta0 = fit$coefficients
  )
  prep <- fastFGEE:::fgee_fastk_prepare(
    st, fx$data, fx$namesd, fit, K = 5L, seed = 8L, exact = TRUE
  )

  d_hold <- vapply(
    prep$folds$sets,
    function(ii) rowSums(st$D_exact[, ii, drop = FALSE]),
    numeric(st$p)
  )
  holdout_rows <- vapply(
    prep$folds$sets,
    function(ii) sum(st$cluster_size[ii]),
    numeric(1)
  )
  expected <- matrix(st$d_exact_sum, st$p, prep$K) - d_hold
  expected <- sweep(
    expected, 2L,
    sum(st$cluster_size) / (sum(st$cluster_size) - holdout_rows),
    "*"
  )
  expect_equal(prep$D, expected, tolerance = 1e-12)

  ev <- fastFGEE:::fgee_fastk_score_grad(
    prep, c(0.8, 1.4), need_gradient = TRUE
  )
  u <- log10(c(0.8, 1.4) / prep$lambda_frem)
  numeric <- finite_diff(function(v) {
    fastFGEE:::fgee_fastk_score_grad(
      prep, prep$lambda_frem * 10^v, need_gradient = FALSE
    )$score
  }, u)
  expect_equal(ev$gradient, numeric, tolerance = 2e-5)
})

test_that("non-Gaussian fastK memory modes give the same score and gradient", {
  fx <- make_working_fixture(N = 12L, n_long = 3L, n_fun = 5L, p = 7L,
                             seed = 215)
  fit <- make_fake_fit(p = 7L, q = 2L, family = binomial())
  fx$data[, Y := rbinom(.N, 1L, 0.45)]
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores"
  )
  lambda <- c(0.6, 1.9)
  ans <- lapply(c("speed", "balanced", "lowmem"), function(mem) {
    prep <- fastFGEE:::fgee_fastk_prepare(
      st, fx$data, fx$namesd, fit, K = 4L, seed = 8L, memory = mem
    )
    fastFGEE:::fgee_fastk_score_grad(prep, lambda, need_gradient = TRUE)
  })
  expect_equal(ans[[2L]]$score, ans[[1L]]$score, tolerance = 1e-11)
  expect_equal(ans[[3L]]$score, ans[[1L]]$score, tolerance = 1e-11)
  expect_equal(ans[[2L]]$gradient, ans[[1L]]$gradient, tolerance = 1e-10)
  expect_equal(ans[[3L]]$gradient, ans[[1L]]$gradient, tolerance = 1e-10)
})
