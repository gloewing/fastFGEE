# Internal numerical linear algebra -------------------------------------------

#' @keywords internal
#' @noRd
.fgee_validate_square_numeric <- function(x, name = "x") {
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) < 1L || nrow(x) != ncol(x)) {
    stop(name, " must be a non-empty numeric square matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

# Finiteness is checked by checked_symmetric_copy() in
# src/fgee_linear_algebra.cpp, which raises the same condition. Doing it in R
# as well cost a full O(p^2) pass plus a p x p logical on every call -- 24 of
# 158 microseconds at p = 60. The base-R fallback has no compiled guard behind
# it, so it validates for itself, the same split this file already uses for the
# symmetric part.
#' @keywords internal
#' @noRd
.fgee_require_finite <- function(x, name = "x") {
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.fgee_symmetric_part <- function(x) 0.5 * (x + t(x))

#' @keywords internal
#' @noRd
.fgee_restore_matrix_dimnames <- function(x, row_names = NULL,
                                          col_names = NULL) {
  if (!is.null(dim(x))) dimnames(x) <- list(row_names, col_names)
  x
}

#' @keywords internal
#' @noRd
fgee_sympd_inverse <- function(x, fallback = c("base", "error")) {
  fallback <- match.arg(fallback)
  x <- .fgee_validate_square_numeric(x)
  dn <- dimnames(x)
  row_names <- if (is.null(dn)) NULL else dn[[1L]]
  col_names <- if (is.null(dn)) NULL else dn[[2L]]
  ans <- tryCatch(fgee_sympd_inverse_cpp(x), error = identity)
  if (!inherits(ans, "error")) {
    return(.fgee_restore_matrix_dimnames(ans, row_names, col_names))
  }
  if (identical(fallback, "error")) {
    stop(conditionMessage(ans), call. = FALSE)
  }

  .fgee_require_finite(x)
  x <- .fgee_symmetric_part(x)
  base_chol <- tryCatch(chol(x), error = identity)
  if (!inherits(base_chol, "error")) {
    out <- chol2inv(base_chol)
    return(.fgee_restore_matrix_dimnames(out, row_names, col_names))
  }

  warning(
    "Compiled and base Cholesky inversions failed; using base::solve(). ",
    "The matrix may not be positive definite.",
    call. = FALSE
  )
  out <- tryCatch(
    solve(x),
    error = function(e_solve) stop(
      "symmetric matrix inversion failed in the compiled Cholesky, ",
      "base Cholesky, and base solve paths: ", conditionMessage(e_solve),
      call. = FALSE
    )
  )
  .fgee_restore_matrix_dimnames(out, row_names, col_names)
}

#' @keywords internal
#' @noRd
fgee_sympd_solve <- function(x, b, fallback = c("base", "error")) {
  fallback <- match.arg(fallback)
  x <- .fgee_validate_square_numeric(x)
  was_vector <- is.null(dim(b))
  b_names <- if (was_vector) names(b) else NULL
  b_dimnames <- if (was_vector) NULL else dimnames(b)
  bmat <- if (was_vector) matrix(as.numeric(b), ncol = 1L) else as.matrix(b)
  if (!is.numeric(bmat) || nrow(bmat) != nrow(x) || any(!is.finite(bmat))) {
    stop("b must be numeric, finite, and have nrow(b) equal to nrow(x).",
         call. = FALSE)
  }
  storage.mode(bmat) <- "double"

  finish <- function(out) {
    if (was_vector) {
      out <- as.numeric(out[, 1L])
      names(out) <- b_names
      return(out)
    }
    .fgee_restore_matrix_dimnames(
      out,
      if (!is.null(b_dimnames)) b_dimnames[[1L]] else NULL,
      if (!is.null(b_dimnames)) b_dimnames[[2L]] else NULL
    )
  }

  ans <- tryCatch(fgee_sympd_solve_cpp(x, bmat), error = identity)
  if (!inherits(ans, "error")) return(finish(ans))
  if (identical(fallback, "error")) {
    stop(conditionMessage(ans), call. = FALSE)
  }

  .fgee_require_finite(x)
  x <- .fgee_symmetric_part(x)
  R <- tryCatch(chol(x), error = function(e) NULL)
  if (!is.null(R)) {
    return(finish(backsolve(R, forwardsolve(t(R), bmat))))
  }

  warning(
    "Compiled and base Cholesky solves failed; using base::solve(). ",
    "The matrix may not be positive definite.",
    call. = FALSE
  )
  finish(solve(x, bmat))
}
