# Beta regression (logit link) in the exact fastK criterion.

beta_fam <- function(phi = 13.227, carry_precision = TRUE) {
  # `precision` is carried deliberately. Unlike negative binomial, whose theta
  # get_family_info() can recover by parsing "Negative Binomial(2.5)", beta has
  # no string-parse fallback (R/family_fns.R:357-381) and defaults phi to 1.0
  # silently. A family object without it therefore makes the legacy path
  # disagree with the fastK path by an eta-free constant -- which is exactly
  # what carry_precision = FALSE is here to demonstrate.
  f <- list(family = sprintf("Beta regression(%s)", format(phi)),
            link = "logit")
  if (carry_precision) f$precision <- phi
  structure(f, class = "family")
}

beta_prep <- function(phi = 13.227, memory = "balanced", seed = 313, N = 12L,
                      n_long = 3L, n_fun = 5L, p = 6L, K = 4L, sets = NULL,
                      mu = 0.45, y_override = NULL, ...) {
  fx <- make_working_fixture(N = N, n_long = n_long, n_fun = n_fun, p = p,
                             seed = seed)
  set.seed(seed + 5L)
  yy <- if (is.null(y_override)) {
    stats::rbeta(nrow(fx$data), shape1 = mu * phi, shape2 = (1 - mu) * phi)
  } else y_override
  fx$data[, Y := yy]
  rec <- structure(
    list(family = "beta", parameter = "precision", value = phi,
         precision = phi, method = "fixed", source = "test"),
    class = c("fgee_nuisance", "list")
  )
  data.table::setattr(fx$data, "nuisance", rec)
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "full", copy_dt = TRUE
  )
  fit <- make_fake_fit(p = p, nsdf = 1L, q = 2L, family = beta_fam(phi))
  list(prep = fastFGEE:::fgee_fastk_prepare(st, fx$data, fx$namesd, fit,
                                            K = K, seed = 11L, sets = sets,
                                            memory = memory, ...),
       working = st, fit = fit, fx = fx)
}

# How many observations sit on the mu clip at a given lambda. The beta kink is
# reachable -- mu is clipped at 1e-6, i.e. eta = +/-13.8 -- unlike the
# effectively unreachable log-link kinks, so the finite-difference test has to
# be run somewhere interior or it fails depending on which side the step lands.
n_clipped <- function(prep, lam) {
  r <- fastFGEE:::fgee_fastk_score_grad(prep, lam, need_gradient = FALSE,
                                        return_fold = TRUE)
  b <- r$beta_folds
  tot <- 0L
  for (k in seq_len(prep$K)) {
    ii <- prep$fold_starts[k]:prep$fold_ends[k]
    eta <- as.numeric(prep$X_eval[ii, , drop = FALSE] %*% b[, k])
    mu <- stats::plogis(eta)
    tot <- tot + sum(mu <= prep$clip_prob | mu >= 1 - prep$clip_prob)
  }
  tot
}

test_that("beta regression reaches the fastK criterion at all", {
  g <- beta_prep()
  expect_identical(g$prep$family_key, "beta")
  expect_equal(g$prep$nuisance_value, 13.227)
  expect_equal(g$prep$clip_y, 1e-6)
  r <- fastFGEE:::fgee_fastk_score_grad(g$prep, c(1, 1), need_gradient = TRUE)
  expect_true(is.finite(r$score))
  expect_true(all(is.finite(r$gradient)))
})

test_that("a missing precision is refused rather than defaulted", {
  fx <- make_working_fixture(N = 8L, n_long = 3L, n_fun = 5L, p = 5L)
  fx$data[, Y := stats::runif(.N, 0.2, 0.8)]
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "scores", copy_dt = TRUE)
  fit <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L, family = beta_fam())
  expect_error(
    fastFGEE:::fgee_fastk_prepare(st, fx$data, fx$namesd, fit, K = 4L,
                                  seed = 1L),
    "requires a positive precision"
  )
})

test_that("responses outside the unit interval are rejected up front", {
  n <- 12L * 3L * 5L
  y_bad <- stats::runif(n, 0.2, 0.8); y_bad[10L] <- 1.4
  expect_error(beta_prep(y_override = y_bad), "responses in \\[0, 1\\]")
  y_neg <- stats::runif(n, 0.2, 0.8); y_neg[3L] <- -0.01
  expect_error(beta_prep(y_override = y_neg), "responses in \\[0, 1\\]")
  # exactly 0 and 1 are allowed: the loss clips them via clip_y
  y_edge <- stats::runif(n, 0.2, 0.8); y_edge[1L] <- 0; y_edge[2L] <- 1
  expect_silent(invisible(beta_prep(y_override = y_edge)))
})

test_that("the analytic beta gradient matches finite differences", {
  for (phi in c(5, 13.227, 40)) {
    g <- beta_prep(phi = phi)
    for (lam in list(c(0.6, 1.4), c(2.5, 0.35))) {
      expect_identical(n_clipped(g$prep, lam), 0L,
                       info = sprintf("phi=%g lambda=%s clipped", phi,
                                      paste(lam, collapse = ",")))
      u <- log10(lam / g$prep$lambda_frem)
      ga <- fastFGEE:::fgee_fastk_score_grad(
        g$prep, g$prep$lambda_frem * 10^u, TRUE)$gradient
      gn <- finite_diff(function(v) {
        fastFGEE:::fgee_fastk_score_grad(
          g$prep, g$prep$lambda_frem * 10^v, FALSE)$score
      }, u)
      expect_equal(ga, gn, tolerance = 2e-5,
                   info = sprintf("phi=%g lambda=%s", phi,
                                  paste(lam, collapse = ",")))
    }
  }
})

test_that("all three memory modes agree for beta", {
  lam <- c(0.9, 1.6)
  ref <- NULL
  for (mem in c("balanced", "speed", "lowmem")) {
    g <- beta_prep(memory = mem)
    r <- fastFGEE:::fgee_fastk_score_grad(g$prep, lam, need_gradient = TRUE)
    if (is.null(ref)) ref <- r else {
      expect_equal(r$score, ref$score, tolerance = 1e-11, info = mem)
      expect_equal(r$gradient, ref$gradient, tolerance = 1e-10, info = mem)
    }
  }
})

test_that("staged beta fastK reproduces the legacy objective", {
  g <- beta_prep(memory = "speed", K = 4L)
  folds <- fastFGEE:::.fgee_make_cluster_folds(g$working$cluster_id, K = 4L,
                                              seed = 11L)
  cv_grid <- list(c(0.1, 1, 10), c(0.1, 1, 10), c(0.5, 1, 2))
  compact <- fastFGEE:::fgee_tune_fastk_staged(g$prep, cv_grid, verbose = FALSE)
  legacy_grid <- fastFGEE:::.fgee_prepare_staged_grids(cv_grid,
                                                       g$prep$lambda_frem)
  legacy <- fastFGEE:::fun.gee1step.cv(
    w = g$working$W,
    d = fastFGEE:::.fgee_working_dlist(g$working),
    grid = legacy_grid, data = g$fx$data, namesd = g$fx$namesd,
    cname_ = "cname_", fit.initial = g$fit, cv = "fastkfold", K = 4L,
    sets = folds$sets, folds.list = replicate(4L, g$fit$coefficients),
    seed = 11L, exact = FALSE, loss = "nll"
  )
  expect_equal(compact$mse, legacy$mse, tolerance = 5e-10)
  expect_equal(unname(compact$lambda.star), unname(legacy$lambda.star),
               tolerance = 1e-12)
})

test_that("theta and precision are recovered from the family string, or warn", {
  # Both families now read a parenthesised value out of the family string, so a
  # hand-built family carrying its parameter only there is handled correctly
  # rather than silently replaced by 1.
  fit_str <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L,
                           family = beta_fam(13.227, carry_precision = FALSE))
  expect_equal(fastFGEE:::get_family_info(fit_str)$dispersion_cpp, 13.227)

  fit_slot <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L,
                            family = beta_fam(13.227, carry_precision = TRUE))
  expect_equal(fastFGEE:::get_family_info(fit_slot)$dispersion_cpp, 13.227)

  fit_nb <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L,
                          family = structure(list(family = "Negative Binomial(2.5)",
                                                  link = "log"),
                                             class = "family"))
  expect_equal(fastFGEE:::get_family_info(fit_nb)$dispersion_cpp, 2.5)

  # When nothing is recoverable the fallback to 1 is announced, because for beta
  # phi = 1 gives the plausible-looking and entirely wrong v = mu(1-mu)/2.
  fit_none <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L,
                            family = structure(list(family = "Beta regression",
                                                    link = "logit"),
                                               class = "family"))
  expect_warning(fi <- fastFGEE:::get_family_info(fit_none),
                 "precision = 1")
  expect_equal(fi$dispersion_cpp, 1.0)

  fit_nb_none <- make_fake_fit(p = 5L, nsdf = 1L, q = 2L,
                               family = structure(list(family = "Negative Binomial",
                                                       link = "log"),
                                                  class = "family"))
  expect_warning(fi2 <- fastFGEE:::get_family_info(fit_nb_none),
                 "theta = 1")
  expect_equal(fi2$dispersion_cpp, 1.0)
})
