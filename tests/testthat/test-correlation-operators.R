test_that("exact exchangeable inverse matches dense solve", {
  set.seed(1)
  for (n in c(2L, 5L, 20L)) {
    for (rho in c(-0.02, 0, 0.25, 0.8)) {
      if (rho <= -1 / (n - 1)) next
      B <- matrix(rnorm(n * 4), n, 4)
      R <- matrix(rho, n, n); diag(R) <- 1
      got <- fastFGEE:::.fgee_apply_exchangeable_inverse(B, rho)
      expect_equal(got, solve(R, B), tolerance = 1e-10)
    }
  }
})

test_that("exact regular AR1 inverse matches dense solve", {
  set.seed(2)
  for (n in c(2L, 5L, 20L)) {
    for (rho in c(-0.7, 0, 0.25, 0.9)) {
      B <- matrix(rnorm(n * 4), n, 4)
      R <- toeplitz(rho^(0:(n - 1L)))
      got <- fastFGEE:::.fgee_apply_ar1_regular_inverse(B, rho)
      expect_equal(got, solve(R, B), tolerance = 2e-10)
    }
  }
})

test_that("exact solvers agree with SuperGauss", {
  skip_if_not_installed("SuperGauss")
  set.seed(3)
  for (corr in c("ar1", "exchangeable")) {
    n <- 13L; rho <- 0.55
    B <- matrix(rnorm(n * 5), n, 5)
    ex <- fastFGEE:::.fgee_apply_corr_inverse(
      B, corr = corr, rho = rho, solver = "exact"
    )
    sg <- fastFGEE:::.fgee_apply_corr_inverse(
      B, corr = corr, rho = rho, solver = "supergauss"
    )
    expect_equal(ex, sg, tolerance = 1e-8)
  }
})

test_that("batched Kronecker inverse matches dense solve", {
  set.seed(4)
  nf <- 7L; nl <- 4L; q <- 5L
  rf <- 0.6; rl <- 0.3
  Q <- matrix(rnorm(nf * nl * q), nf * nl, q)
  af <- function(B) fastFGEE:::.fgee_apply_ar1_regular_inverse(B, rf)
  al <- function(B) fastFGEE:::.fgee_apply_exchangeable_inverse(B, rl)
  got <- fastFGEE:::.fgee_apply_kron_inverse(Q, nf, nl, af, al)
  Rf <- toeplitz(rf^(0:(nf - 1L)))
  Rl <- matrix(rl, nl, nl); diag(Rl) <- 1
  expect_equal(got, solve(kronecker(Rl, Rf), Q), tolerance = 5e-10)
})

test_that("exact correlation operators remain stable near valid boundaries", {
  set.seed(5)

  # Negative exchangeable correlation close to the positive-definite boundary.
  n <- 8L
  rho_ex <- -0.14
  B <- matrix(rnorm(n * 3L), n, 3L)
  R_ex <- matrix(rho_ex, n, n)
  diag(R_ex) <- 1
  expect_equal(
    fastFGEE:::.fgee_apply_exchangeable_inverse(B, rho_ex),
    solve(R_ex, B),
    tolerance = 2e-9
  )

  # AR(1) permits negative rho as long as |rho| < 1.
  rho_ar <- -0.95
  R_ar <- toeplitz(rho_ar^(0:(n - 1L)))
  expect_equal(
    fastFGEE:::.fgee_apply_ar1_regular_inverse(B, rho_ar),
    solve(R_ar, B),
    tolerance = 2e-8
  )
})

test_that("generic SuperGauss ACF hook matches a dense Toeplitz solve", {
  skip_if_not_installed("SuperGauss")
  set.seed(6)
  n <- 17L
  h <- 0:(n - 1L)
  a <- sqrt(3) / 4
  acf <- (1 + a * h) * exp(-a * h) # Matern 3/2 on a regular grid
  B <- matrix(rnorm(n * 4L), n, 4L)

  got <- fastFGEE:::.fgee_apply_supergauss_acf_inverse(B, acf = acf)
  expect_equal(got, solve(toeplitz(acf), B), tolerance = 1e-7)
})

test_that("FPCA Woodbury inverse matches the dense covariance inverse", {
  set.seed(7)
  n <- 12L
  r <- 3L
  Phi <- qr.Q(qr(matrix(rnorm(n * r), n, r)))
  fpca <- list(
    efunctions = Phi,
    evalues = c(2.2, 0.8, 0.3),
    sigma2 = 0.4,
    argvals = seq_len(n)
  )
  B <- matrix(rnorm(n * 5L), n, 5L)
  op <- fastFGEE:::.fgee_make_fpca_inverse_operator(fpca)
  V <- Phi %*% diag(fpca$evalues) %*% t(Phi) + fpca$sigma2 * diag(n)

  expect_equal(op(B, seq_len(n)), solve(V, B), tolerance = 2e-9)
})

test_that("invalid correlation parameters fail before backend fallback", {
  B <- matrix(1, 5L, 2L)
  expect_error(
    fastFGEE:::.fgee_apply_corr_inverse(
      B, corr = "ar1", rho = 1, solver = "auto"
    ),
    "requires \\|rho\\| < 1"
  )
  expect_error(
    fastFGEE:::.fgee_apply_corr_inverse(
      B, corr = "exchangeable", rho = -0.3, solver = "auto"
    ),
    "not positive definite"
  )
})
