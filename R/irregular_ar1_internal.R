# Internal irregular-time AR(1) helpers ---------------------------------------
#
# Independent implementation of the Gaussian Markov precision and profiled
# likelihood in Allevius (2018), "On the precision matrix of an irregularly
# sampled AR(1) process". No source code from the archived R package is used.

#' @keywords internal
#' @noRd
.fgee_order_iar1 <- function(residual, time) {
  residual <- as.numeric(residual)
  time <- as.numeric(time)
  if (length(residual) != length(time)) {
    stop("residual and time must have equal length.", call. = FALSE)
  }
  keep <- is.finite(residual) & is.finite(time)
  residual <- residual[keep]
  time <- time[keep]
  if (length(residual) < 2L) {
    stop("At least two finite irregular AR(1) observations are required.",
         call. = FALSE)
  }
  ord <- order(time)
  residual <- residual[ord]
  time <- time[ord]
  if (anyDuplicated(time)) {
    stop("Irregular AR(1) observation times must be unique.", call. = FALSE)
  }
  if (any(diff(time) <= 0)) {
    stop("Irregular AR(1) observation times must be strictly increasing.",
         call. = FALSE)
  }
  list(residual = residual, time = time)
}

#' @keywords internal
#' @noRd
fgee_iar1_precision_bands <- function(time, rho) {
  time <- as.numeric(time)
  if (any(!is.finite(time))) {
    stop("time must contain only finite values.", call. = FALSE)
  }
  if (length(time) > 1L && any(diff(time) <= 0)) {
    stop("time must be strictly increasing.", call. = FALSE)
  }
  fgee_iar1_precision_bands_cpp(time, as.numeric(rho)[1L])
}

#' @keywords internal
#' @noRd
fgee_iar1_precision <- function(time, rho, sparse = FALSE) {
  time <- as.numeric(time)
  rho <- as.numeric(rho)[1L]
  if (!isTRUE(sparse)) return(fgee_iar1_precision_cpp(time, rho))
  bands <- fgee_iar1_precision_bands_cpp(time, rho)
  if (length(time) == 1L) {
    return(Matrix::Diagonal(n = 1L, x = bands$diagonal))
  }
  Matrix::bandSparse(
    n = length(time),
    k = c(-1L, 0L, 1L),
    diagonals = list(bands$off_diagonal, bands$diagonal,
                     bands$off_diagonal),
    symmetric = FALSE
  )
}

#' @keywords internal
#' @noRd
fgee_iar1_apply_precision <- function(rhs, time, rho) {
  was_vector <- is.null(dim(rhs))
  rhs_names <- if (was_vector) names(rhs) else NULL
  rhs_dimnames <- if (was_vector) NULL else dimnames(rhs)
  rhs <- if (was_vector) matrix(as.numeric(rhs), ncol = 1L) else as.matrix(rhs)
  if (!is.numeric(rhs) || any(!is.finite(rhs))) {
    stop("rhs must contain only finite numeric values.", call. = FALSE)
  }
  storage.mode(rhs) <- "double"
  out <- fgee_iar1_apply_precision_cpp(
    rhs,
    as.numeric(time),
    as.numeric(rho)[1L]
  )
  if (was_vector) {
    out <- as.numeric(out[, 1L])
    names(out) <- rhs_names
    return(out)
  }
  dimnames(out) <- rhs_dimnames
  out
}

#' @keywords internal
#' @noRd
fgee_iar1_profile_nll <- function(rho, residual, time) {
  dat <- .fgee_order_iar1(residual, time)
  fgee_iar1_profile_nll_cpp(
    dat$residual,
    dat$time,
    as.numeric(rho)[1L]
  )
}

# Estimate a common one-unit rho from one series, a matrix of series, or a list
# of series. This is retained as an internal method-development helper; the
# validated public estimator continues to use its established correlation
# update schedule.
#' @keywords internal
#' @noRd
fgee_estimate_iar1 <- function(
    residual,
    time,
    lower = 1e-6,
    upper = 1 - 1e-6,
    weighting = c("cluster", "observation"),
    control = list()) {
  weighting <- match.arg(weighting)
  lower <- as.numeric(lower)[1L]
  upper <- as.numeric(upper)[1L]
  if (!is.finite(lower) || !is.finite(upper) ||
      lower <= 0 || upper >= 1 || lower >= upper) {
    stop("lower and upper must satisfy 0 < lower < upper < 1.", call. = FALSE)
  }

  make_series <- function(residual, time) {
    if (is.list(residual) && !is.data.frame(residual)) {
      if (is.list(time)) {
        if (length(time) != length(residual)) {
          stop("residual and time lists differ in length.", call. = FALSE)
        }
        return(Map(.fgee_order_iar1, residual, time))
      }
      return(lapply(residual, .fgee_order_iar1, time = time))
    }
    if (is.matrix(residual)) {
      if (length(time) == ncol(residual)) {
        return(lapply(seq_len(nrow(residual)), function(i) {
          .fgee_order_iar1(residual[i, ], time)
        }))
      }
      if (length(time) == nrow(residual)) {
        return(lapply(seq_len(ncol(residual)), function(i) {
          .fgee_order_iar1(residual[, i], time)
        }))
      }
      stop("One matrix dimension must equal length(time).", call. = FALSE)
    }
    list(.fgee_order_iar1(residual, time))
  }

  series <- make_series(residual, time)
  if (!length(series)) {
    stop("At least one irregular AR(1) residual series is required.",
         call. = FALSE)
  }
  sizes <- vapply(series, function(z) length(z$residual), integer(1L))

  # Each series criterion is an innovation-variance-profiled negative twice
  # log likelihood and therefore grows with its number of observations.  Work
  # with a per-observation criterion first.  Equal-cluster weighting then gives
  # every series equal influence, while observation weighting reconstructs the
  # pooled average criterion sum(q_i) / sum(n_i).
  weights <- if (identical(weighting, "cluster")) {
    rep(1 / length(series), length(series))
  } else {
    sizes / sum(sizes)
  }

  objective <- function(rho) {
    values <- vapply(series, function(z) {
      fgee_iar1_profile_nll_cpp(z$residual, z$time, rho)
    }, numeric(1L))
    sum(weights * (values / sizes))
  }

  tol <- if (!is.null(control$tol)) {
    as.numeric(control$tol)[1L]
  } else {
    .Machine$double.eps^0.25
  }
  if (!is.finite(tol) || tol <= 0) {
    stop("control$tol must be a positive finite scalar.", call. = FALSE)
  }
  opt <- stats::optimize(objective, interval = c(lower, upper), tol = tol)
  boundary_tol <- max(10 * tol, sqrt(.Machine$double.eps))
  boundary <- (opt$minimum - lower <= boundary_tol) ||
    (upper - opt$minimum <= boundary_tol)

  structure(
    list(
      rho = as.numeric(opt$minimum),
      objective = as.numeric(opt$objective),
      convergence = 0L,
      converged = TRUE,
      boundary = boundary,
      bounds = c(lower = lower, upper = upper),
      weighting = weighting,
      n_series = length(series),
      n_used = sum(sizes),
      message = if (boundary) {
        "rho estimate is at or near an optimization boundary"
      } else {
        NULL
      }
    ),
    class = c("fgee_iar1_estimate", "list")
  )
}
