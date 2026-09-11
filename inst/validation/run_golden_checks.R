suppressPackageStartupMessages({
  library(fastFGEE)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else getwd()
source(file.path(script_dir, "helpers.R"))

cat("Running fastFGEE optimized-engine golden checks\n")

set.seed(1)
for (n in c(2L, 8L, 25L)) {
  B <- matrix(rnorm(n * 4L), n, 4L)
  for (rho in c(-0.7, 0, 0.75)) {
    R <- toeplitz(rho^(0:(n - 1L)))
    stopifnot(isTRUE(all.equal(
      fastFGEE:::.fgee_apply_ar1_regular_inverse(B, rho),
      solve(R, B), tolerance = 2e-8
    )))
  }
  lower <- -1 / (n - 1L)
  for (rho in c(max(lower + 0.02, -0.2), 0, 0.7)) {
    R <- matrix(rho, n, n); diag(R) <- 1
    stopifnot(isTRUE(all.equal(
      fastFGEE:::.fgee_apply_exchangeable_inverse(B, rho),
      solve(R, B), tolerance = 2e-8
    )))
  }
}
cat("  correlation operators: PASS\n")

fx <- validation_fixture(N = 4L, n_long = 5L, n_fun = 7L, p = 5L)
for (cl in c("independent", "ar1", "exchangeable")) {
  for (cf in c("independent", "ar1", "exchangeable")) {
    oldW <- fastFGEE:::.getW(
      fx$data, fx$namesd, "cname_", corr_fn = cf, corr_long = cl,
      ensure_order = "setorderv"
    )
    oldD <- fastFGEE:::.getD(
      fx$data, fx$namesd, "cname_", corr_fn = cf, corr_long = cl,
      resid_col = "resid", ensure_order = "setorderv"
    )
    nw <- fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, corr_fn = cf, corr_long = cl,
      retain = "full", corr_solver = "exact"
    )
    stopifnot(isTRUE(all.equal(
      lapply(nw$W, unname),
      lapply(oldW, unname),
      tolerance = 2e-8
    )))
    stopifnot(isTRUE(all.equal(
      lapply(fastFGEE:::.fgee_working_dlist(nw), unname),
      lapply(oldD, unname),
      tolerance = 2e-8
    )))
  }
}
cat("  one-pass W/D versus legacy: PASS\n")

fit <- validation_fake_fit(p = 7L, q = 2L, family = binomial())
fb <- validation_fixture(N = 15L, n_long = 3L, n_fun = 5L, p = 7L,
                         family = "binomial", seed = 4L)
st <- fastFGEE:::fgee_build_working_stats(fb$data, fb$namesd, retain = "scores")
fp <- fastFGEE:::fgee_fastk_prepare(st, fb$data, fb$namesd, fit, K = 5L, seed = 9L)
u <- log10(c(0.8, 1.6) / fp$lambda_frem)
ga <- fastFGEE:::fgee_fastk_score_grad(
  fp, fp$lambda_frem * 10^u, need_gradient = TRUE
)$gradient
gn <- validation_finite_diff(function(v) {
  fastFGEE:::fgee_fastk_score_grad(
    fp, fp$lambda_frem * 10^v, need_gradient = FALSE
  )$score
}, u)
stopifnot(max(abs(ga - gn)) < 2e-5)
cat("  analytic fastK gradient: PASS\n")

qp <- fastFGEE:::fgee_qreml_prepare(st, fit)
qi <- fastFGEE:::fgee_qreml_information(qp)
phi <- min(max(qi$phi_pen, 1), 8)
ne <- st$N / phi
qa <- fastFGEE:::fgee_qreml_score_grad(
  qp, qp$lambda_frem * 10^u, ne, need_gradient = TRUE
)$gradient
qn <- validation_finite_diff(function(v) {
  fastFGEE:::fgee_qreml_score_grad(
    qp, qp$lambda_frem * 10^v, ne, need_gradient = FALSE
  )$score
}, u)
stopifnot(max(abs(qa - qn)) < 2e-5)
cat("  sandwich qREML gradient: PASS\n")

stopifnot(identical(
  fastFGEE:::.fgee_resolve_working_retain(
    "auto", joint.CI = "wild", var.type = "sandwich"
  ),
  "scores"
))
cat("  default wild-CI retention: PASS\n")
cat("All golden checks passed.\n")
