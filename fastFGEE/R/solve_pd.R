#' @keywords internal
#' @noRd
solve_pd <- function(a, b = NULL) {
  # Try sanic first (fastest for PD matrices)
  if (requireNamespace("sanic", quietly = TRUE)) {
    if (is.null(b)) {
      return(sanic::solve_chol(a = a))
    } else {
      return(sanic::solve_chol(a = a, b = b))
    }
  }

  # Try Cholesky-based methods (much faster than solve for PD matrices)
  chol_A <- tryCatch(chol(a), error = function(e) NULL)

  if (!is.null(chol_A)) {
    if (is.null(b)) {
      # Matrix inversion
      return(chol2inv(chol_A))
    } else {
      # Solve A %*% X = B via Cholesky
      return(backsolve(chol_A, forwardsolve(t(chol_A), b)))
    }
  }

  # Fall back to solve if Cholesky fails
  warning("Using base::solve() - matrix may not be positive definite")
  if (is.null(b)) {
    return(solve(a))
  } else {
    return(solve(a, b))
  }
}
