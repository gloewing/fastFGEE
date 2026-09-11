test_that("registered SPD inverse and solve agree with base R", {
  set.seed(9006)
  Z <- matrix(rnorm(48), 8, 6)
  A <- crossprod(Z) + diag(0.5, 6)
  B <- matrix(rnorm(18), 6, 3)
  dimnames(A) <- list(paste0("r", seq_len(6)), paste0("r", seq_len(6)))
  dimnames(B) <- list(rownames(A), paste0("b", seq_len(3)))

  inv_cpp <- fastFGEE:::fgee_sympd_inverse(A, fallback = "error")
  sol_cpp <- fastFGEE:::fgee_sympd_solve(A, B, fallback = "error")

  expect_equal(inv_cpp, solve(A), tolerance = 1e-10)
  expect_equal(sol_cpp, solve(A, B), tolerance = 1e-10)
  expect_equal(unname(A %*% inv_cpp), diag(6), tolerance = 1e-10)
  expect_equal(A %*% sol_cpp, B, tolerance = 1e-10)
  expect_identical(dimnames(inv_cpp), dimnames(A))
  expect_identical(dimnames(sol_cpp), dimnames(B))

  A_integer <- matrix(c(4L, 1L, 1L, 3L), 2, 2)
  b_integer <- matrix(c(1L, 2L), 2, 1)
  # base::solve() drops dimnames where these routines preserve them by design,
  # so compare values unnamed; preservation is asserted separately above.
  expect_equal(
    unname(fastFGEE:::fgee_sympd_solve(A_integer, b_integer, fallback = "error")),
    unname(solve(A_integer, b_integer)),
    tolerance = 1e-12
  )
})

test_that("irregular AR1 precision agrees with a dense correlation inverse", {
  times <- c(0, 0.2, 0.75, 1.9, 2.1, 4.6)
  rho <- 0.73
  R <- outer(times, times, function(a, b) rho^abs(a - b))
  Q <- fastFGEE:::fgee_iar1_precision(times, rho)
  expect_equal(Q, solve(R), tolerance = 1e-9)

  X <- cbind(seq_along(times), sin(times), rep(1, length(times)))
  expect_equal(
    fastFGEE:::fgee_iar1_apply_precision(X, times, rho),
    solve(R, X),
    tolerance = 1e-9
  )
  expect_equal(
    fastFGEE:::.fgee_apply_corr_inverse(
      X, corr = "ar1", rho = rho, times = times,
      grid_type = "irregular", solver = "exact"
    ),
    solve(R, X),
    tolerance = 1e-9
  )
})

test_that("irregular operator reduces to regular AR1 on unit gaps", {
  times <- 0:7
  rho <- 0.41
  rhs <- matrix(seq_len(24), 8, 3)
  expect_equal(
    fastFGEE:::fgee_iar1_apply_precision(rhs, times, rho),
    fastFGEE:::.fgee_apply_ar1_regular_inverse(rhs, rho),
    tolerance = 1e-11
  )
})

test_that("profile objective agrees with the dense Gaussian expression", {
  residual <- c(-0.7, 0.1, 1.2, -0.2, 0.5, -1.1)
  times <- c(0, 0.4, 1.1, 1.8, 3.2, 4.9)
  for (rho in c(0.15, 0.55, 0.9)) {
    R <- outer(times, times, function(a, b) rho^abs(a - b))
    q <- drop(crossprod(residual, solve(R, residual)))
    expected <- length(residual) * log(q / length(residual)) +
      as.numeric(determinant(R, logarithm = TRUE)$modulus)
    observed <- fastFGEE:::fgee_iar1_profile_nll(rho, residual, times)
    expect_equal(observed, expected, tolerance = 1e-9)
  }
})

test_that("irregular AR1 input validation is explicit", {
  expect_error(
    fastFGEE:::fgee_iar1_precision(c(0, 1, 1), 0.5),
    "strictly increasing"
  )
  expect_error(
    fastFGEE:::fgee_iar1_precision(c(0, 1, 2), -0.1),
    "0 <= rho < 1"
  )
  expect_error(
    fastFGEE:::fgee_estimate_iar1(list(), list()),
    "At least one"
  )
  expect_error(
    fastFGEE:::fgee_estimate_iar1(
      c(0, 1), c(0, 1), control = list(tol = -1)
    ),
    "positive finite"
  )
})

test_that("single-time irregular precision is the scalar identity", {
  dense <- fastFGEE:::fgee_iar1_precision(2.5, 0.7)
  sparse <- fastFGEE:::fgee_iar1_precision(2.5, 0.7, sparse = TRUE)
  expect_equal(dense, matrix(1, 1, 1), tolerance = 0)
  expect_equal(as.matrix(sparse), matrix(1, 1, 1), tolerance = 0)
})

test_that("rho zero is the irregular working-independence boundary", {
  times <- c(0, 0.1, 0.8, 2.5)
  rhs <- matrix(seq_len(8), 4, 2)
  expect_equal(
    fastFGEE:::fgee_iar1_precision(times, 0),
    diag(length(times)),
    tolerance = 0
  )
  expect_equal(
    fastFGEE:::fgee_iar1_apply_precision(rhs, times, 0),
    rhs,
    tolerance = 0
  )
  expect_equal(
    fastFGEE:::fgee_iar1_apply_precision(matrix(1:8, 4, 2), times, 0),
    matrix(as.numeric(1:8), 4, 2),
    tolerance = 0
  )
})


test_that("irregular AR1 common-rho weighting has an explicit estimand", {
  residual <- list(
    c(-0.7, 0.2, 1.1, -0.3),
    c(0.1, -0.4, 0.8, 0.5, -0.2, 0.3, -0.9)
  )
  times <- list(
    c(0, 0.4, 1.7, 3.1),
    c(0, 0.2, 0.9, 1.4, 2.8, 4.0, 5.5)
  )
  sizes <- vapply(residual, length, integer(1L))
  qvals <- function(rho) vapply(seq_along(residual), function(i) {
    fastFGEE:::fgee_iar1_profile_nll(rho, residual[[i]], times[[i]])
  }, numeric(1L))
  cluster_objective <- function(rho) mean(qvals(rho) / sizes)
  observation_objective <- function(rho) sum(qvals(rho)) / sum(sizes)

  ref_cluster <- optimize(cluster_objective, c(1e-6, 1 - 1e-6))
  ref_observation <- optimize(observation_objective, c(1e-6, 1 - 1e-6))
  fit_cluster <- fastFGEE:::fgee_estimate_iar1(
    residual, times, weighting = "cluster"
  )
  fit_observation <- fastFGEE:::fgee_estimate_iar1(
    residual, times, weighting = "observation"
  )

  expect_equal(fit_cluster$rho, ref_cluster$minimum, tolerance = 1e-8)
  expect_equal(fit_cluster$objective, ref_cluster$objective, tolerance = 1e-8)
  expect_equal(fit_observation$rho, ref_observation$minimum, tolerance = 1e-8)
  expect_equal(
    fit_observation$objective,
    ref_observation$objective,
    tolerance = 1e-8
  )
  expect_identical(fit_cluster$n_series, 2L)
  expect_identical(fit_cluster$n_used, sum(sizes))
})

test_that("SPD helpers fail explicitly or use the documented base fallback", {
  indefinite <- diag(c(1, -2))
  rhs <- c(a = 3, b = -4)

  expect_error(
    fastFGEE:::fgee_sympd_inverse(indefinite, fallback = "error"),
    "positive definite"
  )
  expect_error(
    fastFGEE:::fgee_sympd_solve(indefinite, rhs, fallback = "error"),
    "positive definite"
  )

  expect_warning(
    inverse <- fastFGEE:::fgee_sympd_inverse(indefinite, fallback = "base"),
    "base::solve"
  )
  expect_warning(
    solution <- fastFGEE:::fgee_sympd_solve(
      indefinite, rhs, fallback = "base"
    ),
    "base::solve"
  )
  expect_equal(unname(inverse), unname(solve(indefinite)), tolerance = 1e-12)
  expect_equal(unname(solution), unname(solve(indefinite, rhs)),
               tolerance = 1e-12)
  expect_identical(names(solution), names(rhs))
})
