# Optional rank-one Gram fast path --------------------------------------------
#
# .fgee_build_working_stats_core() consumes the working inverse only through
# crossprod(Z, R^{-1} Q) -- the m x k matrix R^{-1} Q is never needed in its own
# right.  When exactly one axis carries an exchangeable correlation and the
# other is independent, the Gram matrix Q' R^{-1} Q can therefore be assembled
# from rank-one sufficient statistics, without ever applying R^{-1}.
#
# The exchangeable inverse on an axis of length n is a * I - b * J with
# a = 1/(1-rho) and b = rho/((1-rho)(1+(n-1)rho)).  Two cases follow:
#
#   A. corr_fn = independent, corr_long = exchangeable.  rho_long may differ
#      from one functional point to the next -- which is the usual case, since
#      it is estimated per functional point -- so R^{-1} is block diagonal with
#      its own (a_f, b_f) per block:
#        Q'R^{-1}Q = Q' diag(a_f(i)) Q - S' diag(b) S,  S = per-functional sums
#
#   B. corr_long = independent, corr_fn = exchangeable.  Mirror image, with
#      rho_fn allowed to differ across longitudinal observations:
#        Q'R^{-1}Q = Q' diag(a_l(i)) Q - S' diag(b) S,  S = per-longitudinal sums
#
# a is strictly positive (1 - rho > tol), so the first term is formed as a
# symmetric rank-k update of the scaled design -- crossprod(Q * sqrt(a)) -- which
# is half the arithmetic of the general crossprod(Q, Q * a) and is exactly
# symmetric by construction.
#
# Case A in particular replaces a loop of n_fun separate correlation solves per
# cluster (.fgee_apply_cluster_inverse()'s split-by-functional-index branch)
# with two BLAS calls.
#
# Deliberately NOT handled here:
#
#   * AR(1) and FPCA axes -- their inverses have no identity-plus-low-rank form,
#     and reducing only the other axis measured no faster than the existing
#     kernel.
#   * exchangeable x exchangeable -- the compiled two-axis inverse kernel
#     already handles it, and the Gram route measured neutral end to end
#     (0.97x), so it is left where it was.
#
# This is an optimisation only: the operator is identical to the one
# .fgee_apply_cluster_inverse() applies, positive-definiteness is checked with
# the same rule and tolerance as .fgee_validate_corr_parameter(), and every case
# that does not qualify returns NULL so the caller falls through to the
# reference path.
#
# Disable with options(fastFGEE.corr.gram = FALSE).

#' @keywords internal
#' @noRd
fgee_corr_gram_ok <- function() {
  isTRUE(getOption("fastFGEE.corr.gram", TRUE)) &&
    is.function(tryCatch(get(".fgee_gram_axis_sums",
                             envir = asNamespace("fastFGEE")),
                         error = function(e) NULL))
}

# One rho per level of the requested axis, or NULL when rho is not constant
# within a level (in which case no exchangeable block structure exists).
#' @keywords internal
#' @noRd
.fgee_gram_level_rho <- function(rho, n_fun, n_long, by, tol = 1e-12) {
  n <- n_fun * n_long
  r <- suppressWarnings(as.numeric(rho))
  if (length(r) == 1L) r <- rep(r, n)
  if (length(r) != n || !all(is.finite(r))) return(NULL)
  M <- matrix(r, nrow = n_fun, ncol = n_long)
  if (identical(by, "fn")) {
    v <- M[, 1L]
    if (n_long > 1L && max(abs(M - v)) > tol) return(NULL)
  } else {
    v <- M[1L, ]
    if (n_fun > 1L && max(abs(M - rep(v, each = n_fun))) > tol) return(NULL)
  }
  v
}

# Exchangeable inverse written as a * I - b * J, vectorised over rho.
# rho = 0 gives a = 1, b = 0; an axis of length one is the identity.  NULL
# signals "not positive definite", matching .fgee_validate_corr_parameter().
# Note a > 0 always, which is what licenses the sqrt() in the caller.
#' @keywords internal
#' @noRd
.fgee_gram_exch_coef <- function(rho, n, tol = 1e-10) {
  if (n <= 1L) return(list(a = rep(1, length(rho)), b = rep(0, length(rho))))
  if (!all(is.finite(rho))) return(NULL)
  d1 <- 1 - rho
  d2 <- 1 + (n - 1) * rho
  if (any(d1 <= tol) || any(d2 <= tol)) return(NULL)
  list(a = 1 / d1, b = rho / (d1 * d2))
}

# Returns Q' R^{-1} Q, or NULL to fall through to the apply path.
#' @keywords internal
#' @noRd
.fgee_cluster_gram <- function(Q, idx_fn, idx_long, corr_fn, corr_long,
                               rho_fn = NULL, rho_long = NULL,
                               grid_type = "auto", corr_solver = "auto",
                               index_fn = "yindex.vec", index_long = "time") {
  if (!fgee_corr_gram_ok()) return(NULL)
  # An explicit request for the SuperGauss backend is honoured, not bypassed.
  if (identical(corr_solver, "supergauss")) return(NULL)
  if (!is.matrix(Q) || !is.double(Q)) return(NULL)

  corr_fn <- .fgee_normalize_corr(corr_fn, allow_fpca = TRUE)
  corr_long <- .fgee_normalize_corr(corr_long, allow_fpca = FALSE)
  # Exactly one exchangeable axis, the other independent.
  exch_fn <- corr_fn == "exchangeable" && corr_long == "independent"
  exch_long <- corr_long == "exchangeable" && corr_fn == "independent"
  if (!exch_fn && !exch_long) return(NULL)

  grid <- tryCatch(
    .fgee_tensor_grid_info(idx_fn, idx_long, index_fn = index_fn,
                           index_long = index_long, require_complete = FALSE),
    error = function(e) NULL
  )
  # An incomplete grid falls through, so the apply path raises the same error
  # it always did when the Kronecker route requires completeness.
  if (is.null(grid) || !isTRUE(grid$complete)) return(NULL)
  n_fun <- grid$n_fun
  n_long <- grid$n_long
  if (nrow(Q) != n_fun * n_long) return(NULL)

  # The compiled sums assume the canonical row order (functional index varying
  # fastest), exactly as the existing Kronecker kernel does.  Verifying it is
  # O(m) and negligible beside the crossproducts.
  key <- match(idx_fn, grid$fn_levels) +
    (match(idx_long, grid$long_levels) - 1L) * n_fun
  if (!identical(key, seq_len(nrow(Q)))) return(NULL)

  if (exch_long) {
    v <- .fgee_gram_level_rho(rho_long, n_fun, n_long, by = "fn")
    if (is.null(v)) return(NULL)
    co <- .fgee_gram_exch_coef(v, n_long)
    if (is.null(co)) return(NULL)
    wt <- rep(sqrt(co$a), times = n_long)
  } else {
    v <- .fgee_gram_level_rho(rho_fn, n_fun, n_long, by = "long")
    if (is.null(v)) return(NULL)
    co <- .fgee_gram_exch_coef(v, n_fun)
    if (is.null(co)) return(NULL)
    wt <- rep(sqrt(co$a), each = n_fun)
  }

  S <- tryCatch(
    .fgee_gram_axis_sums(Q, as.integer(n_fun), as.integer(n_long), exch_long),
    error = function(e) NULL
  )
  if (is.null(S)) return(NULL)

  W <- crossprod(Q * wt) - crossprod(S, S * co$b)

  if (!is.matrix(W) || any(dim(W) != c(ncol(Q), ncol(Q))) || !all(is.finite(W))) {
    return(NULL)
  }
  W
}
