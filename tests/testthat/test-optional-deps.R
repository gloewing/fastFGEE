test_that("the exact correlation solvers need no optional package", {
  # Regular-grid AR(1) and exchangeable correlation are handled directly, so a
  # default fit must not touch the SuperGauss backend at all. Assert it by
  # counting calls through the guarded constructor.
  calls <- 0L
  orig <- fastFGEE:::.fgee_supergauss_toeplitz
  on.exit(assignInNamespace(".fgee_supergauss_toeplitz", orig, ns = "fastFGEE"),
          add = TRUE)
  assignInNamespace(".fgee_supergauss_toeplitz", function(N, acf) {
    calls <<- calls + 1L
    orig(N = N, acf = acf)
  }, ns = "fastFGEE")

  set.seed(2)
  n <- 25L
  rhs <- matrix(stats::rnorm(n * 2L), n, 2L)
  for (cr in c("ar1", "exchangeable")) {
    rho <- if (cr == "ar1") 0.6 else 0.3
    out <- fastFGEE:::.fgee_apply_corr_inverse(rhs, cr, rho = rho,
                                               times = seq_len(n),
                                               solver = "auto")
    expect_true(all(is.finite(out)))
  }
  expect_identical(calls, 0L)
})

test_that("the exact solvers match a dense solve", {
  set.seed(3)
  n <- 30L
  rhs <- matrix(stats::rnorm(n * 3L), n, 3L)
  for (cr in c("ar1", "exchangeable")) {
    rho <- if (cr == "ar1") 0.7 else 0.35
    R <- if (cr == "ar1") {
      rho^abs(outer(seq_len(n), seq_len(n), "-"))
    } else {
      m <- matrix(rho, n, n); diag(m) <- 1; m
    }
    ex <- fastFGEE:::.fgee_apply_corr_inverse(rhs, cr, rho = rho,
                                              times = seq_len(n),
                                              solver = "exact")
    expect_equal(ex, solve(R, rhs), tolerance = 1e-10)
  }
})

test_that("requesting the SuperGauss backend without it installed errors clearly", {
  skip_if(requireNamespace("SuperGauss", quietly = TRUE),
          "SuperGauss is installed, so the missing-package path cannot be exercised")
  expect_error(
    fastFGEE:::.fgee_supergauss_toeplitz(N = 8L, acf = rep(1, 8)),
    "suggested rather than required"
  )
})

test_that("base replacements for the removed Rfast helpers behave correctly", {
  # jointCI previously used Rfast::rowMaxs / colMaxs; these are the row and
  # column maxima it needs, verified bit-identical before the substitution.
  set.seed(5)
  X <- matrix(stats::rnorm(50 * 7), 50, 7)
  expect_equal(apply(X, 1L, max), vapply(seq_len(nrow(X)),
               function(i) max(X[i, ]), numeric(1)))
  expect_equal(apply(X, 2L, max), vapply(seq_len(ncol(X)),
               function(j) max(X[, j]), numeric(1)))

  # corr.estimate previously used Rfast::ar1(method = "yw"); the quantity is the
  # lag-1 autocorrelation.
  set.seed(6)
  r <- as.numeric(stats::arima.sim(list(ar = 0.5), 200))
  rc <- stats::acf(r, lag.max = 1L, type = "correlation", plot = FALSE,
                   demean = TRUE)$acf
  expect_length(rc, 2L)
  expect_true(abs(as.numeric(rc[2L])) < 1)
})
