test_that("compact sandwich equals list-based sandwich", {
  set.seed(121)
  p <- 6L; n <- 20L
  W <- lapply(seq_len(n), function(i) crossprod(matrix(rnorm(p * p), p, p)) + diag(p))
  D <- matrix(rnorm(p * n), p, n)
  P <- diag(c(0, rep(0.2, p - 1L)))
  st <- structure(list(
    N = n, p = p, W_bar = Reduce(`+`, W) / n, W_sum = Reduce(`+`, W),
    D = D, W = W, D_exact = NULL, d_bar = rowMeans(D),
    dd_sum = tcrossprod(D), cluster_id = as.character(seq_len(n))
  ), class = c("fgee_working_stats", "list"))
  old <- fastFGEE:::var.est(
    di = lapply(seq_len(n), function(i) D[, i]), wi = W,
    beta2 = rep(0, p), penalty_diag = P, var.type = "sandwich"
  )
  new <- fastFGEE:::fgee_var_from_stats(
    working = st, beta = rep(0, p), penalty_diag = P,
    var.type = "sandwich"
  )
  expect_equal(new, old, tolerance = 1e-12)
})

test_that("compact exact-Gaussian sandwich equals the list formulation", {
  set.seed(123)
  p <- 6L; n <- 18L
  beta <- rnorm(p)
  W <- lapply(seq_len(n), function(i) {
    crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  })
  U <- matrix(rnorm(p * n), p, n)
  D_exact <- do.call(cbind, lapply(seq_len(n), function(i) {
    U[, i] + as.numeric(W[[i]] %*% beta)
  }))
  P <- diag(c(0, rep(0.15, p - 1L)))
  st <- structure(list(
    N = n,
    p = p,
    W_bar = Reduce(`+`, W) / n,
    W_sum = Reduce(`+`, W),
    D = U,
    W = NULL,
    D_exact = D_exact,
    d_bar = rowMeans(U),
    dd_sum = tcrossprod(U),
    cluster_id = as.character(seq_len(n)),
    retain = "scores"
  ), class = c("fgee_working_stats", "list"))

  old <- fastFGEE:::var.est(
    di = lapply(seq_len(n), function(i) D_exact[, i]),
    wi = W,
    beta2 = beta,
    penalty_diag = P,
    var.type = "sandwich",
    exact = TRUE
  )
  new <- fastFGEE:::fgee_var_from_stats(
    working = st,
    beta = beta,
    penalty_diag = P,
    var.type = "sandwich",
    exact = TRUE
  )

  expect_equal(new, old, tolerance = 1e-12)
})

test_that("compact wild bootstrap matches list-based implementation", {
  set.seed(122)
  p <- 5L; n <- 16L
  W <- lapply(seq_len(n), function(i) crossprod(matrix(rnorm(p * p), p, p)) + diag(p))
  D <- matrix(rnorm(p * n), p, n)
  P <- diag(c(0, rep(0.3, p - 1L)))
  beta <- rnorm(p)
  A <- list(a = matrix(rnorm(25), 5, p), b = matrix(rnorm(20), 4, p))
  # Equalize row counts, as required by the package implementation.
  A$b <- rbind(A$b, A$b[1L, ])
  Sigma <- diag(p) / 10

  old <- fastFGEE:::wild_studentized_boot_crit(
    beta_hat = beta, Sigma_hat = Sigma,
    di = lapply(seq_len(n), function(i) D[, i]), wi = W,
    penalty_diag = P, A_list = A, B = 100L, seed = 9L,
    progress_every = 0L, beta_center = beta
  )
  new <- fastFGEE:::wild_studentized_boot_crit_compact(
    beta_hat = beta, Sigma_hat = Sigma, D = D,
    Wbar = Reduce(`+`, W) / n, penalty_diag = P,
    A_list = A, B = 100L, seed = 9L, progress_every = 0L,
    beta_center = beta
  )
  expect_equal(new$joint, old$joint, tolerance = 1e-12)
  expect_equal(new$pointwise_crit, old$pointwise_crit, tolerance = 1e-12)
  expect_equal(new$ci$df, old$ci$df, tolerance = 1e-12)
})


test_that("aggregate retention is sufficient for sandwich covariance", {
  set.seed(124)
  p <- 5L; n <- 22L
  W <- lapply(seq_len(n), function(i) {
    crossprod(matrix(rnorm(p * p), p, p)) + diag(p)
  })
  D <- matrix(rnorm(p * n), p, n)
  P <- diag(c(0, rep(0.25, p - 1L)))
  st <- structure(list(
    N = n,
    p = p,
    W_bar = Reduce(`+`, W) / n,
    W_sum = Reduce(`+`, W),
    D = NULL,
    W = NULL,
    d_bar = rowMeans(D),
    dd_sum = tcrossprod(D),
    cluster_id = as.character(seq_len(n)),
    retain = "aggregate"
  ), class = c("fgee_working_stats", "list"))

  old <- fastFGEE:::var.est(
    di = lapply(seq_len(n), function(i) D[, i]),
    wi = W,
    beta2 = rep(0, p),
    penalty_diag = P,
    var.type = "sandwich"
  )
  new <- fastFGEE:::fgee_var_from_stats(
    working = st,
    beta = rep(0, p),
    penalty_diag = P,
    var.type = "sandwich"
  )
  expect_equal(new, old, tolerance = 2e-12)
})

test_that("precomputed eta and fitted values update the embedded model without MM", {
  mod <- list(
    beta = c(0.2, -0.1),
    rho = list(),
    model = list(
      coefficients = c(0, 0),
      family = gaussian(),
      linear.predictors = numeric(3),
      fitted.values = numeric(3)
    )
  )
  eta <- c(-1, 0, 2)
  mu <- eta
  got <- fastFGEE:::fgee_model_update(mod, eta = eta, fitted = mu)
  expect_equal(got$model$coefficients, mod$beta)
  expect_equal(got$model$linear.predictors, eta)
  expect_equal(got$model$fitted.values, mu)
})
