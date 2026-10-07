# gee.fit = FALSE: cluster sandwich around the initial pffr fit ---------------
#
# mgcv's smoothing parameters penalise the total (unscaled) deviance or RSS,
# while the one-step bread is on the averaged scale Wbar + P.  These tests
# compare the gee.fit = FALSE covariance against a cluster sandwich built
# directly from the mgcv fit (bread from fit$Vp, scores from fit$y and
# fit$fitted.values), using the package's conventions (centred cluster scores,
# factor N / (N - 1)).

sim_initial_fit_data <- function(family, N = 30L, n_i = 3L, L = 20L, seed = 1L) {
  set.seed(seed)
  s <- seq(0, 1, length.out = L)
  n <- N * n_i
  ID <- rep(seq_len(N), each = n_i)
  X1 <- stats::rnorm(n)
  u <- matrix(stats::rnorm(2L * N), N, 2L)
  eta <- outer(rep(1, n), 0.5 + 0.5 * sin(2 * pi * s)) +
    X1 %o% (0.5 * cos(2 * pi * s)) +
    0.5 * u[ID, 1L] %o% rep(1, L) + 0.5 * u[ID, 2L] %o% s
  Y <- switch(family,
    gaussian = eta + matrix(stats::rnorm(n * L, sd = 0.8), n, L),
    poisson = matrix(stats::rpois(n * L, exp(eta)), n, L)
  )
  d <- data.frame(ID = ID, time = rep(seq_len(n_i), N), X1 = X1)
  d$Y <- I(Y)
  list(data = d, s = s)
}

fit_initial_pffr <- function(sim, family) {
  d <- sim$data
  s <- sim$s
  refund::pffr(
    Y ~ X1, yind = s, data = d, family = family,
    algorithm = "bam", discrete = TRUE, sandwich = "none",
    bs.yindex = list(bs = "ps", k = 8, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 8, m = c(2, 1))
  )
}

fit_initial_sandwich <- function(sim, fit, family, ...) {
  suppressMessages(fastFGEE:::.fgee_fit_internal(
    Y ~ X1, data = sim$data, cluster = "ID", family = family, time = "time",
    pffr.mod = fit, gee.fit = FALSE, var.type = "sandwich",
    bs = "ps", knots = 7, ...
  ))
}

# Reference cluster sandwich around the mgcv fit.  For the Gaussian family the
# package's working variance v(s) is the residual variance at each grid point,
# so the reference uses the same weights: bread X'V^-1 X + S_lambda / sig2.
reference_sandwich <- function(sim, fit) {
  X <- suppressWarnings(stats::predict(fit, type = "lpmatrix"))
  L <- length(sim$s)
  cl <- sim$data$ID[rep(seq_len(nrow(sim$data)), each = L)]
  res <- fit$y - fit$fitted.values
  if (identical(fit$family$family, "gaussian")) {
    S_lambda <- fit$sig2 * solve(fit$Vp) - crossprod(X)
    grid_id <- rep(seq_len(L), nrow(sim$data))
    w <- 1 / stats::ave(res, grid_id, FUN = stats::var)
    H <- crossprod(X * w, X) + S_lambda / fit$sig2
  } else {
    # Canonical link: Vp = (X'WX + S_lambda)^-1, score X'(y - mu).
    w <- 1
    H <- solve(fit$Vp)
  }
  scores <- rowsum(X * (res * w), cl)
  scores <- sweep(scores, 2L, colMeans(scores))
  N <- nrow(scores)
  Hinv <- solve(H)
  N / (N - 1) * Hinv %*% crossprod(scores) %*% Hinv
}

se_ratio <- function(V, V_ref) sqrt(diag(V)) / sqrt(diag(V_ref))

test_that("gee.fit = FALSE sandwich matches the cluster sandwich around the pffr fit", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")

  for (fam_name in c("gaussian", "poisson")) {
    family <- switch(fam_name, gaussian = stats::gaussian(), poisson = stats::poisson())
    sim <- sim_initial_fit_data(fam_name)
    fit <- fit_initial_pffr(sim, family)
    V_ref <- reference_sandwich(sim, fit)
    for (engine in c("optimized", "legacy")) {
      ans <- fit_initial_sandwich(
        sim, fit, family, joint.CI = FALSE, working.engine = engine,
        sp.method = if (engine == "legacy") "legacy" else "auto"
      )
      expect_equal(unname(ans$beta), unname(fit$coefficients), tolerance = 1e-10)
      expect_equal(unname(ans$vb), unname(V_ref), tolerance = 1e-6,
                   label = paste(fam_name, engine, "vb"))

      # Returned smoothing parameters are on the averaged scale.
      scale <- if (fam_name == "gaussian") 30 * fit$sig2 else 30
      expect_equal(as.numeric(ans$lambda), as.numeric(fit$sp) / scale,
                   tolerance = 1e-10)
    }
  }
})

test_that("gee.fit = FALSE penalty solves the pffr estimating equation", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")

  # At the penalised fit, mean cluster score = P beta on the averaged scale
  # (exact for a canonical link and known scale).
  sim <- sim_initial_fit_data("poisson", seed = 2L)
  fit <- fit_initial_pffr(sim, stats::poisson())
  ans <- fit_initial_sandwich(sim, fit, stats::poisson(), joint.CI = FALSE)
  dbar <- ans$working0$d_bar
  expect_lt(
    sqrt(sum((dbar - ans$pen.mat %*% ans$beta)^2)) / sqrt(sum(dbar^2)),
    1e-6
  )

  # The EDF used for the wild-bootstrap t adjustment reproduces mgcv's edf.
  edf <- fastFGEE:::fgee_effective_df_from_Wbar(
    Wbar = ans$working0$W_bar,
    penalty_mat = ans$pen.mat,
    A_list = list(diag(length(ans$beta)))
  )
  expect_equal(unname(edf$diagS), unname(fit$edf), tolerance = 1e-5)
})

test_that("gee.fit = FALSE wild-bootstrap CIs use the averaged-scale penalty", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")

  sim <- sim_initial_fit_data("poisson", seed = 3L)
  fit <- fit_initial_pffr(sim, stats::poisson())
  ans <- fit_initial_sandwich(sim, fit, stats::poisson(), joint.CI = "wild")
  info <- ans$model$crit$info

  # Term EDFs equal mgcv's edf summed over each term's coefficients.
  by <- vapply(fit$smooth, function(sm) sm$by, character(1))
  cols <- function(sm) sm$first.para:sm$last.para
  edf_mgcv <- c(
    sum(fit$edf[c(1L, cols(fit$smooth[[which(by == "NA")]]))]),
    sum(fit$edf[cols(fit$smooth[[which(by == "X1")]])])
  )
  expect_equal(as.numeric(info$edf_r), edf_mgcv, tolerance = 1e-5)
  expect_equal(as.numeric(info$df_r), 30 - edf_mgcv, tolerance = 1e-5)

  # The intervals are studentized with the corrected sandwich covariance.
  expect_equal(unname(ans$model$Vp), unname(reference_sandwich(sim, fit)),
               tolerance = 1e-6)
  expect_true(all(is.finite(ans$model$ci$df$X1$se)))
})

test_that("gee.fit = FALSE standard errors do not shrink relative to the reference as N grows", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")

  ratios <- vapply(c(20L, 120L), function(N) {
    sim <- sim_initial_fit_data("poisson", N = N, seed = 4L)
    fit <- fit_initial_pffr(sim, stats::poisson())
    ans <- fit_initial_sandwich(sim, fit, stats::poisson(), joint.CI = FALSE)
    stats::median(se_ratio(ans$vb, reference_sandwich(sim, fit)))
  }, numeric(1))
  expect_equal(ratios, c(1, 1), tolerance = 1e-6)
})
