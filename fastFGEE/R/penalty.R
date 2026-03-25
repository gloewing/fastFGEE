#' @keywords internal
#' @noRd
penalty_setup <- function(model, unpenalized = NULL) {
  if (is.null(unpenalized)) unpenalized <- model$nsdf
  unpenalized <- as.integer(unpenalized)

  p <- length(model$coefficients)
  if (!is.finite(p) || p <= 0) stop("Could not determine coefficient length from model$coefficients.")

  comps <- list()
  ncomp_by_smooth <- integer(length(model$smooth))

  for (i in seq_along(model$smooth)) {
    sm <- model$smooth[[i]]
    ind <- sm$first.para:sm$last.para
    if (length(ind) == 0) next

    ncomp <- length(sm$S)
    ncomp_by_smooth[i] <- ncomp
    if (ncomp == 0) next

    for (k in seq_len(ncomp)) {
      Sblk <- as.matrix(sm$S[[k]])
      if (!all(dim(Sblk) == length(ind))) {
        stop(sprintf(
          "Penalty dim mismatch for smooth %d comp %d: dim(S)=%s, #ind=%d",
          i, k, paste(dim(Sblk), collapse = "x"), length(ind)
        ))
      }
      comps[[length(comps) + 1L]] <- list(ind = ind, S = Sblk, smooth = i, comp = k)
    }
  }

  list(
    p = p,
    unpenalized = unpenalized,
    comps = comps,
    n_sm = length(model$smooth),
    n_pen = length(comps),
    ncomp_by_smooth = ncomp_by_smooth
  )
}

#' @keywords internal
#' @noRd
penalty_from_setup <- function(setup, lambda) {
  p <- setup$p
  nsdf <- setup$unpenalized
  comps <- setup$comps

  n_sm  <- setup$n_sm
  n_pen <- setup$n_pen
  ncomp_by_smooth <- setup$ncomp_by_smooth

  lambda <- as.numeric(lambda)

  if (length(lambda) == 1L) {
    lambda_pen <- rep(lambda, n_pen)
  } else if (length(lambda) == n_sm) {
    lambda_pen <- rep(lambda, times = ncomp_by_smooth)
  } else if (length(lambda) == n_pen) {
    lambda_pen <- lambda
  } else {
    stop(sprintf(
      "lambda must have length 1, #smooths=%d, or #penalties=%d. Got %d.",
      n_sm, n_pen, length(lambda)
    ))
  }

  P <- matrix(0, p, p)
  if (n_pen > 0) {
    for (j in seq_len(n_pen)) {
      ind <- comps[[j]]$ind
      P[ind, ind] <- P[ind, ind] + lambda_pen[j] * comps[[j]]$S
    }
  }

  if (!is.null(nsdf) && nsdf > 0L) {
    P[1:nsdf, ] <- 0
    P[, 1:nsdf] <- 0
  }

  P
}
