#' @keywords internal
#' @noRd
solve_pd <- function(a, b = NULL) {
  A <- as.matrix(a)
  R <- tryCatch(chol(A), error = function(e) NULL)

  if (!is.null(R)) {
    if (is.null(b)) {
      out <- chol2inv(R)
      dimnames(out) <- dimnames(A)
      return(out)
    }
    out <- backsolve(R, forwardsolve(t(R), b))
    if (is.null(dim(b))) {
      out <- as.numeric(out)
      names(out) <- names(b)
    } else {
      dimnames(out) <- dimnames(b)
    }
    return(out)
  }

  # chol() failed, so the matrix is not numerically positive definite. Say so
  # rather than returning a general-inverse result that looks like an SPD one.
  warning("Cholesky factorization failed in solve_pd(); using base::solve(). ",
          "The matrix may not be positive definite.", call. = FALSE)
  if (is.null(b)) solve(A) else solve(A, b)
}

# Reusable positive-definite factorization ------------------------------------
# Used by smoothing optimizers that need multiple solves and/or a log
# determinant from one Cholesky factorization.
#' @keywords internal
#' @noRd
.fgee_factor_pd <- function(a,
                            ridge = TRUE,
                            ridge_multiplier = 1,
                            max_ridge_tries = 6L) {
  A <- as.matrix(a)
  if (nrow(A) != ncol(A)) stop("a must be square.")
  if (any(!is.finite(A))) stop("a contains non-finite values.")
  A <- 0.5 * (A + t(A))

  R <- tryCatch(chol(A), error = function(e) NULL)
  added <- 0
  if (is.null(R) && isTRUE(ridge)) {
    scale <- max(abs(diag(A)), 1)
    base <- sqrt(.Machine$double.eps) * scale * ridge_multiplier
    for (k in seq_len(as.integer(max_ridge_tries))) {
      added <- base * 10^(k - 1L)
      R <- tryCatch(chol(A + diag(added, nrow(A))), error = function(e) NULL)
      if (!is.null(R)) break
    }
  }
  if (is.null(R)) stop("Matrix is not positive definite after ridge stabilization.")

  solve_fun <- function(b) {
    backsolve(R, forwardsolve(t(R), b))
  }

  list(
    R = R,
    solve = solve_fun,
    logdet = 2 * sum(log(diag(R))),
    ridge = added
  )
}
