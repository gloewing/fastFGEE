# The compiled Kronecker inverse must reproduce the R operator it replaces.
# It is an arithmetic optimisation, not a new estimator, so the whole contract
# is exhaustively testable: agreement over the operator's parameter space, plus
# fall-through wherever the kernel cannot represent the requested operator.

skip_if_no_corr_kernel <- function() {
  testthat::skip_if_not(fastFGEE:::fgee_corr_kernel_ok(),
                        "compiled correlation kernel unavailable")
}

kron_reference <- function(Q, n_fun, n_long, corr_fn, rho_fn,
                           corr_long, rho_long) {
  ai <- fastFGEE:::.fgee_apply_corr_inverse
  af <- function(B) ai(B, corr = corr_fn, rho = rho_fn,
                       times = seq_len(n_fun), grid_type = "auto",
                       solver = "exact")
  al <- function(B) ai(B, corr = corr_long, rho = rho_long,
                       times = seq_len(n_long), grid_type = "auto",
                       solver = "exact")
  fastFGEE:::.fgee_apply_kron_inverse(Q, n_fun, n_long, af, al)
}

test_that("the compiled Kronecker inverse matches the R operator", {
  skip_if_no_corr_kernel()
  code <- c(independent = 0L, exchangeable = 1L, ar1 = 2L)
  grid <- expand.grid(
    n_fun = c(1L, 2L, 3L, 8L, 17L),
    n_long = c(1L, 2L, 5L, 12L),
    corr_fn = c("exchangeable", "ar1", "independent"),
    corr_long = c("exchangeable", "ar1", "independent"),
    rho = c(-0.45, 0.05, 0.6, 0.92),
    q = c(1L, 4L),
    stringsAsFactors = FALSE
  )
  worst <- 0
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    m <- g$n_fun * g$n_long
    # An exchangeable axis with negative rho need not be positive definite.
    pd_ok <- function(corr, n) {
      if (corr != "exchangeable" || n <= 1L) return(TRUE)
      1 + (n - 1L) * g$rho > 1e-8
    }
    if (!pd_ok(g$corr_fn, g$n_fun) || !pd_ok(g$corr_long, g$n_long)) next
    set.seed(1000L + r)
    Q <- matrix(stats::rnorm(m * g$q), m, g$q)
    ref <- kron_reference(Q, g$n_fun, g$n_long, g$corr_fn, g$rho,
                          g$corr_long, g$rho)
    got <- fastFGEE:::.fgee_kron_kernel_try(
      Q, g$n_fun, g$n_long, g$corr_fn, g$rho, g$corr_long, g$rho,
      times_fn = seq_len(g$n_fun), times_long = seq_len(g$n_long)
    )
    expect_false(is.null(got), label = paste("kernel declined row", r))
    worst <- max(worst, max(abs(ref - got)) / max(1, max(abs(ref))))
  }
  expect_lt(worst, 1e-12)
})

test_that("the kernel declines operators it cannot represent exactly", {
  skip_if_no_corr_kernel()
  Q <- matrix(stats::rnorm(24), 6L, 4L)
  # irregular AR(1) grid
  expect_null(fastFGEE:::.fgee_kron_kernel_try(
    Q, 3L, 2L, "ar1", 0.5, "exchangeable", 0.3,
    times_fn = c(1, 2, 9), times_long = c(1, 2)))
  # FPCA functional correlation
  expect_null(fastFGEE:::.fgee_kron_kernel_try(
    Q, 3L, 2L, "fpca", NULL, "exchangeable", 0.3,
    times_fn = 1:3, times_long = 1:2))
  # non-positive-definite exchangeable
  expect_null(fastFGEE:::.fgee_kron_kernel_try(
    Q, 3L, 2L, "exchangeable", -0.9, "independent", 0,
    times_fn = 1:3, times_long = 1:2))
  # |rho| >= 1 for AR(1)
  expect_null(fastFGEE:::.fgee_kron_kernel_try(
    Q, 3L, 2L, "ar1", 1, "independent", 0,
    times_fn = 1:3, times_long = 1:2))
  # an explicit SuperGauss request is honoured, not bypassed
  expect_null(fastFGEE:::.fgee_kron_kernel_try(
    Q, 3L, 2L, "ar1", 0.5, "exchangeable", 0.3,
    times_fn = 1:3, times_long = 1:2, solver = "supergauss"))
  # and the option disables it
  old <- options(fastFGEE.corr.kernel = FALSE)
  on.exit(options(old), add = TRUE)
  expect_false(fastFGEE:::fgee_corr_kernel_ok())
  expect_null(fastFGEE:::.fgee_kron_kernel_try(
    Q, 3L, 2L, "ar1", 0.5, "exchangeable", 0.3,
    times_fn = 1:3, times_long = 1:2))
})

test_that("working statistics agree with the kernel disabled", {
  skip_if_no_corr_kernel()
  for (sd in c(481L, 17L, 22L)) {
    fx <- make_working_fixture(N = 9L, n_long = 4L, n_fun = 7L, p = 6L,
                               rho_long = 0.4, rho_fun = 0.55, seed = sd)
    for (cc in list(c("ar1", "exchangeable"), c("exchangeable", "ar1"),
                    c("ar1", "independent"), c("independent", "exchangeable"),
                    c("exchangeable", "exchangeable"), c("ar1", "ar1"))) {
      for (ret in c("scores", "full", "aggregate")) {
        build <- function() fastFGEE:::fgee_build_working_stats(
          fx$data, fx$namesd, corr_fn = cc[1L], corr_long = cc[2L],
          retain = ret, copy_dt = TRUE)
        on_ker <- build()
        old <- options(fastFGEE.corr.kernel = FALSE)
        off_ker <- build()
        options(old)
        lab <- paste(sd, cc[1L], cc[2L], ret)
        expect_equal(on_ker$W_sum, off_ker$W_sum, tolerance = 1e-12, label = lab)
        expect_equal(on_ker$d_sum, off_ker$d_sum, tolerance = 1e-12, label = lab)
        expect_equal(on_ker$dd_sum, off_ker$dd_sum, tolerance = 1e-12, label = lab)
        expect_equal(on_ker$D, off_ker$D, tolerance = 1e-12, label = lab)
        expect_equal(on_ker$W, off_ker$W, tolerance = 1e-12, label = lab)
      }
    }
  }
})
