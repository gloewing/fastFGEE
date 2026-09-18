# The rank-one Gram fast path must reproduce the reference operator exactly
# (to floating-point summation order) and must decline every structure it
# cannot represent, so that .fgee_apply_cluster_inverse() stays in charge.

gram_ref <- function(Q, idx_fn, idx_long, corr_fn, corr_long, rho_fn, rho_long) {
  crossprod(Q, fastFGEE:::.fgee_apply_cluster_inverse(
    Q, idx_fn = idx_fn, idx_long = idx_long, corr_fn = corr_fn,
    corr_long = corr_long, rho_fn = rho_fn, rho_long = rho_long))
}
mk_idx <- function(n_fun, n_long) {
  list(fn = rep(seq_len(n_fun), times = n_long),
       lo = rep(seq_len(n_long), each = n_fun))
}
na_of <- function(Q) rep(NA_real_, nrow(Q))

test_that("the Gram fast path matches the apply path for a scalar rho", {
  for (dims in list(c(12, 5, 4), c(7, 3, 6), c(20, 2, 3))) {
    n_fun <- dims[1]; n_long <- dims[2]; k <- dims[3]
    ii <- mk_idx(n_fun, n_long)
    set.seed(11)
    Q <- matrix(rnorm(n_fun * n_long * k), n_fun * n_long, k)

    got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "independent",
             "exchangeable", na_of(Q), rep(0.7, nrow(Q)))
    expect_false(is.null(got))
    expect_equal(got, gram_ref(Q, ii$fn, ii$lo, "independent", "exchangeable",
                               na_of(Q), rep(0.7, nrow(Q))), tolerance = 1e-12)

    got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "exchangeable",
             "independent", rep(0.4, nrow(Q)), na_of(Q))
    expect_false(is.null(got))
    expect_equal(got, gram_ref(Q, ii$fn, ii$lo, "exchangeable", "independent",
                               rep(0.4, nrow(Q)), na_of(Q)), tolerance = 1e-12)
  }
})

test_that("the Gram fast path matches when rho varies across the other axis", {
  n_fun <- 15; n_long <- 6; k <- 5
  ii <- mk_idx(n_fun, n_long)
  set.seed(12)
  Q <- matrix(rnorm(n_fun * n_long * k), n_fun * n_long, k)

  # rho_long estimated separately at each functional point
  rl <- rep(seq(0.05, 0.85, length.out = n_fun), times = n_long)
  got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "independent",
                                       "exchangeable", na_of(Q), rl)
  expect_false(is.null(got))
  expect_equal(got, gram_ref(Q, ii$fn, ii$lo, "independent", "exchangeable",
                             na_of(Q), rl), tolerance = 1e-12)

  # rho_fn estimated separately at each longitudinal observation
  rf <- rep(seq(0.1, 0.8, length.out = n_long), each = n_fun)
  got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "exchangeable",
                                       "independent", rf, na_of(Q))
  expect_false(is.null(got))
  expect_equal(got, gram_ref(Q, ii$fn, ii$lo, "exchangeable", "independent",
                             rf, na_of(Q)), tolerance = 1e-12)
})

test_that("the Gram fast path matches for negative rho up to the bound", {
  # Exchangeable stays positive definite while 1 + (n - 1) rho > 0, so the
  # admissible negative range is rho > -1/(n - 1) on an axis of length n.
  n_fun <- 15; n_long <- 6; k <- 5
  ii <- mk_idx(n_fun, n_long)
  set.seed(21)
  Q <- matrix(rnorm(n_fun * n_long * k), n_fun * n_long, k)

  # longitudinal axis has length 6, so rho > -0.2
  for (r in c(-0.19, -0.1, -1e-8)) {
    rl <- rep(r, nrow(Q))
    got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "independent",
                                         "exchangeable", na_of(Q), rl)
    expect_false(is.null(got), info = paste("rho_long =", r))
    expect_equal(got, gram_ref(Q, ii$fn, ii$lo, "independent", "exchangeable",
                               na_of(Q), rl), tolerance = 1e-12,
                 info = paste("rho_long =", r))
  }

  # functional axis has length 15, so rho > -1/14
  for (r in c(-0.07, -0.03)) {
    rf <- rep(r, nrow(Q))
    got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "exchangeable",
                                         "independent", rf, na_of(Q))
    expect_false(is.null(got), info = paste("rho_fn =", r))
    expect_equal(got, gram_ref(Q, ii$fn, ii$lo, "exchangeable", "independent",
                               rf, na_of(Q)), tolerance = 1e-12,
                 info = paste("rho_fn =", r))
  }

  # and with rho varying across the sign change from one level to the next
  rl <- rep(seq(-0.19, 0.9, length.out = n_fun), times = n_long)
  got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "independent",
                                       "exchangeable", na_of(Q), rl)
  expect_false(is.null(got))
  expect_equal(got, gram_ref(Q, ii$fn, ii$lo, "independent", "exchangeable",
                             na_of(Q), rl), tolerance = 1e-12)
})

test_that("the Gram fast path declines structures it cannot represent", {
  n_fun <- 10; n_long <- 4; k <- 3
  ii <- mk_idx(n_fun, n_long)
  set.seed(13)
  Q <- matrix(rnorm(n_fun * n_long * k), n_fun * n_long, k)
  r <- rep(0.5, nrow(Q))
  decline <- function(...) expect_null(fastFGEE:::.fgee_cluster_gram(...))

  decline(Q, ii$fn, ii$lo, "ar1", "exchangeable", r, r)
  decline(Q, ii$fn, ii$lo, "exchangeable", "ar1", r, r)
  decline(Q, ii$fn, ii$lo, "ar1", "ar1", r, r)
  decline(Q, ii$fn, ii$lo, "fpca", "exchangeable", r, r)
  decline(Q, ii$fn, ii$lo, "independent", "independent", r, r)
  # both axes exchangeable is left to the compiled two-axis inverse kernel
  decline(Q, ii$fn, ii$lo, "exchangeable", "exchangeable", r, r)
  # rho not constant within a level: no exchangeable block structure exists
  decline(Q, ii$fn, ii$lo, "independent", "exchangeable",
          na_of(Q), runif(nrow(Q), 0.2, 0.8))
  # not positive definite at one level, from either side of the bound
  decline(Q, ii$fn, ii$lo, "independent", "exchangeable", na_of(Q),
          rep(rep(c(0.5, 1), length.out = n_fun), times = n_long))
  decline(Q, ii$fn, ii$lo, "independent", "exchangeable", na_of(Q),
          rep(-1 / (n_long - 1), nrow(Q)))
  decline(Q, ii$fn, ii$lo, "independent", "exchangeable", na_of(Q),
          rep(-0.9, nrow(Q)))
  decline(Q, ii$fn, ii$lo, "independent", "exchangeable", na_of(Q), na_of(Q))
  # an incomplete tensor grid leaves the Kronecker route to the apply path
  decline(Q[-1, , drop = FALSE], ii$fn[-1], ii$lo[-1], "independent",
          "exchangeable", na_of(Q)[-1], r[-1])
  # rows not in the canonical (functional index fastest) order
  perm <- c(2:nrow(Q), 1L)
  decline(Q, ii$fn[perm], ii$lo[perm], "independent", "exchangeable",
          na_of(Q), r)
})

test_that("the fast path can be switched off and leaves results unchanged", {
  n_fun <- 9; n_long <- 4; k <- 3
  ii <- mk_idx(n_fun, n_long)
  set.seed(14)
  Q <- matrix(rnorm(n_fun * n_long * k), n_fun * n_long, k)
  rl <- rep(0.6, nrow(Q))
  on_ <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "independent",
                                       "exchangeable", na_of(Q), rl)
  opt <- options(fastFGEE.corr.gram = FALSE)
  expect_null(fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "independent",
                                            "exchangeable", na_of(Q), rl))
  options(opt)
  expect_equal(on_, gram_ref(Q, ii$fn, ii$lo, "independent", "exchangeable",
                             na_of(Q), rl), tolerance = 1e-12)
})

test_that("the weighted term is formed as an exactly symmetric rank-k update", {
  # a = 1/(1 - rho) is strictly positive, so crossprod(Q * sqrt(a)) is valid
  # and symmetric to the last bit; only the low-rank correction is general.
  n_fun <- 11; n_long <- 5; k <- 4
  ii <- mk_idx(n_fun, n_long)
  set.seed(15)
  Q <- matrix(rnorm(n_fun * n_long * k), n_fun * n_long, k)
  got <- fastFGEE:::.fgee_cluster_gram(Q, ii$fn, ii$lo, "independent",
           "exchangeable", na_of(Q), rep(0.55, nrow(Q)))
  expect_false(is.null(got))
  expect_lt(max(abs(got - t(got))) / max(abs(got)), 1e-14)
})

test_that("the axis-sum kernel forms only the requested axis", {
  n_fun <- 6; n_long <- 4; k <- 3
  set.seed(16)
  Q <- matrix(rnorm(n_fun * n_long * k), n_fun * n_long, k)
  by_fn <- fastFGEE:::.fgee_gram_axis_sums(Q, n_fun, n_long, TRUE)
  by_long <- fastFGEE:::.fgee_gram_axis_sums(Q, n_fun, n_long, FALSE)
  expect_equal(dim(by_fn), c(n_fun, k))
  expect_equal(dim(by_long), c(n_long, k))
  A <- array(Q, dim = c(n_fun, n_long, k))
  expect_equal(by_fn, apply(A, c(1, 3), sum))
  expect_equal(by_long, apply(A, c(2, 3), sum))
})
