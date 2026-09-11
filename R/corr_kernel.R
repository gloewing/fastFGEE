# Optional compiled Kronecker inverse -----------------------------------------
#
# The compiled kernel in src/corr_kernel.cpp is an optimisation of
# .fgee_apply_kron_inverse() only: same operator, same result to machine
# precision, one allocation instead of a dozen.  It is reached only when every
# guard below holds, and .fgee_apply_kron_inverse() remains the reference
# implementation and the fallback -- so a package built without compiled code,
# or with the kernel disabled, produces the same numbers by the same route it
# always did.
#
# Disable with options(fastFGEE.corr.kernel = FALSE).

#' @keywords internal
#' @noRd
fgee_corr_kernel_ok <- function() {
  isTRUE(getOption("fastFGEE.corr.kernel", TRUE)) &&
    is.function(tryCatch(get(".fgee_kron_inverse_kernel",
                             envir = asNamespace("fastFGEE")),
                         error = function(e) NULL))
}

# Correlation code for the kernel, or NA when the kernel cannot represent the
# requested operator (FPCA, irregular AR(1), an out-of-range rho).
#' @keywords internal
#' @noRd
.fgee_corr_kernel_code <- function(corr, rho, n, times, grid_type, tol = 1e-10) {
  if (is.null(corr) || length(corr) != 1L) return(NA_integer_)
  if (identical(corr, "independent") || n <= 1L) return(0L)
  if (!identical(corr, "exchangeable") && !identical(corr, "ar1")) {
    return(NA_integer_)
  }
  rho <- suppressWarnings(as.numeric(rho)[1L])
  if (length(rho) != 1L || !is.finite(rho)) return(NA_integer_)
  if (identical(corr, "exchangeable")) {
    if (1 - rho <= tol || 1 + (n - 1L) * rho <= tol) return(NA_integer_)
    return(1L)
  }
  # AR(1): the closed-form tridiagonal precision is the regular-grid operator.
  if (!.fgee_is_regular_grid(times, grid_type = grid_type)) return(NA_integer_)
  if (1 - rho^2 <= tol) return(NA_integer_)
  2L
}

# Returns the kernel result, or NULL to fall through to the R path.
#' @keywords internal
#' @noRd
.fgee_kron_kernel_try <- function(Q, n_fun, n_long,
                                  corr_fn, rho_fn, corr_long, rho_long,
                                  times_fn, times_long,
                                  grid_type = "auto", solver = "auto",
                                  tol = 1e-10) {
  if (!fgee_corr_kernel_ok()) return(NULL)
  # An explicit request for the SuperGauss backend is honoured, not bypassed.
  if (identical(solver, "supergauss")) return(NULL)
  if (!is.matrix(Q) || !is.double(Q)) return(NULL)
  if (nrow(Q) != as.numeric(n_fun) * as.numeric(n_long)) return(NULL)

  code_fn <- .fgee_corr_kernel_code(corr_fn, rho_fn, n_fun, times_fn,
                                    grid_type, tol)
  if (is.na(code_fn)) return(NULL)
  code_long <- .fgee_corr_kernel_code(corr_long, rho_long, n_long, times_long,
                                      grid_type, tol)
  if (is.na(code_long)) return(NULL)
  if (code_fn == 0L && code_long == 0L) return(Q)

  rf <- if (code_fn == 0L) 0 else as.numeric(rho_fn)[1L]
  rl <- if (code_long == 0L) 0 else as.numeric(rho_long)[1L]
  out <- tryCatch(
    .fgee_kron_inverse_kernel(Q, as.integer(n_fun), as.integer(n_long),
                              as.integer(code_fn), rf,
                              as.integer(code_long), rl),
    error = function(e) NULL
  )
  if (is.null(out) || !is.matrix(out) || any(dim(out) != dim(Q))) return(NULL)
  out
}
