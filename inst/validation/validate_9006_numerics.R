#!/usr/bin/env Rscript

# Focused validation for the numerical and public-interface changes introduced
# in fastFGEE 0.3.0.9006. Run after installing the candidate package.

stopifnot(requireNamespace("fastFGEE", quietly = TRUE))
stopifnot(requireNamespace("Matrix", quietly = TRUE))

cat("fastFGEE version:", as.character(utils::packageVersion("fastFGEE")), "\n")
stopifnot(utils::packageVersion("fastFGEE") == "0.3.0.9006")

# Public estimator surface ----------------------------------------------------
fml <- names(formals(fastFGEE::fgee))
removed <- c("exact", "gee.fit", "max.iter", "tune.method", "working.engine")
stopifnot(!any(removed %in% fml))
stopifnot(identical(
  eval(formals(fastFGEE::fgee)$sp.method),
  c("auto", "fastk_staged", "fastk_grad", "fastk_grad_fast",
    "sandwich_qreml", "qreml_fastk")
))
stopifnot(identical(
  fastFGEE:::.fgee_resolve_sp_method("auto", stats::gaussian(), "identity"),
  "fastk_staged"
))
stopifnot(identical(
  fastFGEE:::.fgee_resolve_sp_method("auto", stats::binomial(), "logit"),
  "fastk_grad_fast"
))

# Registered routines ---------------------------------------------------------
dll <- getLoadedDLLs()[["fastFGEE"]]
stopifnot(!is.null(dll))
registered <- names(getDLLRegisteredRoutines(dll)$.Call)
expected <- c(
  "_fastFGEE_fgee_kron_inverse_kernel",
  "_fastFGEE_fastk_fold_kernel",
  "_fastFGEE_fgee_sympd_inverse_cpp",
  "_fastFGEE_fgee_sympd_solve_cpp",
  "_fastFGEE_fgee_iar1_precision_bands_cpp",
  "_fastFGEE_fgee_iar1_precision_cpp",
  "_fastFGEE_fgee_iar1_apply_precision_cpp",
  "_fastFGEE_fgee_iar1_profile_nll_cpp"
)
stopifnot(setequal(registered, expected))

# Symmetric positive-definite inverse and solve -------------------------------
set.seed(9006)
Z <- matrix(rnorm(48), 8, 6)
A <- crossprod(Z) + diag(0.5, 6)
B <- matrix(rnorm(18), 6, 3)
Ai <- fastFGEE:::fgee_sympd_inverse(A, fallback = "error")
AB <- fastFGEE:::fgee_sympd_solve(A, B, fallback = "error")
stopifnot(max(abs(Ai - solve(A))) < 1e-10)
stopifnot(max(abs(AB - solve(A, B))) < 1e-10)

# Exact irregular-time AR(1) precision ---------------------------------------
times <- c(0, 0.2, 0.75, 1.9, 2.1, 4.6)
rho <- 0.73
R <- outer(times, times, function(a, b) rho^abs(a - b))
Q <- fastFGEE:::fgee_iar1_precision(times, rho)
stopifnot(max(abs(Q - solve(R))) < 1e-9)
X <- cbind(seq_along(times), sin(times), rep(1, length(times)))
stopifnot(max(abs(
  fastFGEE:::fgee_iar1_apply_precision(X, times, rho) - solve(R, X)
)) < 1e-9)
stopifnot(max(abs(
  fastFGEE:::.fgee_apply_corr_inverse(
    X, corr = "ar1", rho = rho, times = times,
    grid_type = "irregular", solver = "exact"
  ) - solve(R, X)
)) < 1e-9)
stopifnot(identical(
  unname(fastFGEE:::fgee_iar1_precision(times, 0)),
  unname(diag(length(times)))
))

# Profile objective agrees with the dense Gaussian expression.
e <- c(-0.7, 0.1, 1.2, -0.2, 0.5, -1.1)
for (rr in c(0.15, 0.55, 0.9)) {
  RR <- outer(times, times, function(a, b) rr^abs(a - b))
  quad <- drop(crossprod(e, solve(RR, e)))
  expected_obj <- length(e) * log(quad / length(e)) +
    as.numeric(determinant(RR, logarithm = TRUE)$modulus)
  observed_obj <- fastFGEE:::fgee_iar1_profile_nll(rr, e, times)
  stopifnot(abs(observed_obj - expected_obj) < 1e-9)
}



# Equal-series versus observation-weighted common-rho criteria have explicit
# interpretations.  This helper remains internal and is not used by the public
# correlation update schedule in this development candidate.
resid_list <- list(
  c(-0.7, 0.2, 1.1, -0.3),
  c(0.1, -0.4, 0.8, 0.5, -0.2, 0.3, -0.9)
)
time_list <- list(
  c(0, 0.4, 1.7, 3.1),
  c(0, 0.2, 0.9, 1.4, 2.8, 4.0, 5.5)
)
sizes <- vapply(resid_list, length, integer(1L))
qvals <- function(rr) vapply(seq_along(resid_list), function(i) {
  fastFGEE:::fgee_iar1_profile_nll(rr, resid_list[[i]], time_list[[i]])
}, numeric(1L))
ref_cluster <- optimize(function(rr) mean(qvals(rr) / sizes),
                        c(1e-6, 1 - 1e-6))
ref_observation <- optimize(function(rr) sum(qvals(rr)) / sum(sizes),
                            c(1e-6, 1 - 1e-6))
fit_cluster <- fastFGEE:::fgee_estimate_iar1(
  resid_list, time_list, weighting = "cluster"
)
fit_observation <- fastFGEE:::fgee_estimate_iar1(
  resid_list, time_list, weighting = "observation"
)
stopifnot(abs(fit_cluster$rho - ref_cluster$minimum) < 1e-8)
stopifnot(abs(fit_cluster$objective - ref_cluster$objective) < 1e-8)
stopifnot(abs(fit_observation$rho - ref_observation$minimum) < 1e-8)
stopifnot(abs(fit_observation$objective - ref_observation$objective) < 1e-8)

cat("Focused 0.3.0.9006 numerical validation passed.\n")
