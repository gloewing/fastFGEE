# Negative binomial (NB2, log link) in the exact fastK criterion.

nb_fam <- function(theta = 2.5) {
  # The literal fitted-style string, not mgcv::nb(): an unfitted nb() reports
  # the clean "negative binomial" and its getTheta(TRUE) silently returns 1, so
  # a constructor-built fixture would exercise neither the parsing nor a real
  # theta.
  structure(list(family = sprintf("Negative Binomial(%s)", format(theta)),
                 link = "log"), class = "family")
}

nb_prep <- function(theta = 2.5, memory = "balanced", seed = 271, N = 12L,
                    n_long = 3L, n_fun = 5L, p = 6L, K = 4L, sets = NULL,
                    mu_scale = 0) {
  fx <- make_working_fixture(N = N, n_long = n_long, n_fun = n_fun, p = p,
                             seed = seed)
  set.seed(seed + 7L)
  mu <- exp(mu_scale)
  fx$data[, Y := stats::rnbinom(.N, size = theta, mu = mu)]
  rec <- structure(
    list(family = "negbinomial", parameter = "theta", value = theta,
         theta = theta, method = "fixed", source = "test"),
    class = c("fgee_nuisance", "list")
  )
  data.table::setattr(fx$data, "nuisance", rec)
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "full", copy_dt = TRUE
  )
  fit <- make_fake_fit(p = p, nsdf = 1L, q = 2L, family = nb_fam(theta))
  list(prep = fastFGEE:::fgee_fastk_prepare(st, fx$data, fx$namesd, fit,
                                            K = K, seed = 11L, sets = sets,
                                            memory = memory),
       working = st, fit = fit, fx = fx)
}

test_that("negative binomial reaches the fastK criterion at all", {
  g <- nb_prep()
  expect_identical(g$prep$family_key, "negbinomial")
  expect_equal(g$prep$nuisance_value, 2.5)
  expect_identical(g$prep$nuisance_source, "working$nuisance")
  r <- fastFGEE:::fgee_fastk_score_grad(g$prep, c(1, 1), need_gradient = TRUE)
  expect_true(is.finite(r$score))
  expect_true(all(is.finite(r$gradient)))
})

test_that("a missing theta is refused rather than defaulted", {
  fx <- make_working_fixture(N = 8L, n_long = 3L, n_fun = 5L, p = 5L)
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "scores", copy_dt = TRUE)
  fit <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L, family = nb_fam(2.5))
  expect_error(
    fastFGEE:::fgee_fastk_prepare(st, fx$data, fx$namesd, fit, K = 4L,
                                  seed = 1L),
    "requires a positive theta"
  )
})

test_that("the analytic negative-binomial gradient matches finite differences", {
  # The returned gradient is with respect to u = log10(lambda / lambda_fREML),
  # so the finite difference has to be taken in that same space.
  for (theta in c(0.8, 2.5, 40)) {
    g <- nb_prep(theta = theta)
    for (lam in list(c(0.4, 1.7), c(3.1, 0.2))) {
      u <- log10(lam / g$prep$lambda_frem)
      ga <- fastFGEE:::fgee_fastk_score_grad(
        g$prep, g$prep$lambda_frem * 10^u, TRUE)$gradient
      gn <- finite_diff(function(v) {
        fastFGEE:::fgee_fastk_score_grad(
          g$prep, g$prep$lambda_frem * 10^v, FALSE)$score
      }, u)
      expect_equal(ga, gn, tolerance = 2e-5,
                   info = sprintf("theta=%g lambda=%s", theta,
                                  paste(lam, collapse = ",")))
    }
  }
})

test_that("the gradient tends to the Poisson gradient as theta grows", {
  # theta (mu - y)/(mu + theta) -> mu - y.  The score carries an extra eta-free
  # (theta + y) log(theta) term that swamps everything at large theta, so the
  # limit is asserted on the gradient, where that term does not appear.
  lam <- c(1.3, 0.6)
  pois <- local({
    fx <- make_working_fixture(N = 12L, n_long = 3L, n_fun = 5L, p = 6L,
                               seed = 271)
    set.seed(278L)
    fx$data[, Y := stats::rnbinom(.N, size = 1e8, mu = 1)]
    st <- fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
      retain = "full", copy_dt = TRUE)
    fit <- make_fake_fit(p = 6L, nsdf = 1L, q = 2L, family = stats::poisson())
    prep <- fastFGEE:::fgee_fastk_prepare(st, fx$data, fx$namesd, fit, K = 4L,
                                          seed = 11L, memory = "balanced")
    fastFGEE:::fgee_fastk_score_grad(prep, lam, need_gradient = TRUE)$gradient
  })
  prev <- Inf
  for (theta in c(1e4, 1e6, 1e8)) {
    g <- nb_prep(theta = theta, seed = 271)
    gn <- fastFGEE:::fgee_fastk_score_grad(g$prep, lam,
                                           need_gradient = TRUE)$gradient
    d <- max(abs(gn - pois)) / max(1, max(abs(pois)))
    expect_lt(d, prev)          # monotone approach
    prev <- d
  }
  expect_lt(prev, 1e-6)
})

test_that("all three memory modes agree for negative binomial", {
  lam <- c(0.9, 2.2)
  ref <- NULL
  for (mem in c("balanced", "speed", "lowmem")) {
    g <- nb_prep(memory = mem)
    r <- fastFGEE:::fgee_fastk_score_grad(g$prep, lam, need_gradient = TRUE)
    if (is.null(ref)) {
      ref <- r
    } else {
      expect_equal(r$score, ref$score, tolerance = 1e-11, info = mem)
      expect_equal(r$gradient, ref$gradient, tolerance = 1e-10, info = mem)
    }
  }
})

test_that("staged negative-binomial fastK reproduces the legacy objective", {
  # An independent implementation of the same loss, written by a different
  # mechanism (R/cv.R:346). The fixture keeps mu in a range where neither the
  # 1e-12 fastK floor nor cv.R's 1e-10 floor binds, so the differing clip
  # constants are provably inert here.
  theta <- 2.5
  g <- nb_prep(theta = theta, memory = "speed", K = 4L)
  folds <- fastFGEE:::.fgee_make_cluster_folds(g$working$cluster_id, K = 4L,
                                              seed = 11L)
  cv_grid <- list(c(0.1, 1, 10), c(0.1, 1, 10), c(0.5, 1, 2))
  compact <- fastFGEE:::fgee_tune_fastk_staged(g$prep, cv_grid, verbose = FALSE)
  legacy_grid <- fastFGEE:::.fgee_prepare_staged_grids(cv_grid,
                                                       g$prep$lambda_frem)
  legacy <- fastFGEE:::fun.gee1step.cv(
    w = g$working$W,
    d = fastFGEE:::.fgee_working_dlist(g$working),
    grid = legacy_grid,
    data = g$fx$data,
    namesd = g$fx$namesd,
    cname_ = "cname_",
    fit.initial = g$fit,
    cv = "fastkfold",
    K = 4L,
    sets = folds$sets,
    folds.list = replicate(4L, g$fit$coefficients),
    seed = 11L,
    exact = FALSE,
    loss = "nll"
  )
  expect_equal(compact$mse, legacy$mse, tolerance = 5e-10)
  expect_equal(unname(compact$lambda.star), unname(legacy$lambda.star),
               tolerance = 1e-12)
})
