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

# Penalty component matrices aligned with the supplied lambda vector ----------
#' @keywords internal
#' @noRd
penalty_components_from_setup <- function(setup, lambda_length) {
  q <- as.integer(lambda_length)[1L]
  if (!is.finite(q) || q < 1L) stop("lambda_length must be a positive integer.")

  p <- setup$p
  n_sm <- setup$n_sm
  n_pen <- setup$n_pen
  comps <- setup$comps

  make_full <- function(ind, S) {
    out <- matrix(0, p, p)
    out[ind, ind] <- S
    if (!is.null(setup$unpenalized) && setup$unpenalized > 0L) {
      out[seq_len(setup$unpenalized), ] <- 0
      out[, seq_len(setup$unpenalized)] <- 0
    }
    out
  }

  if (q == 1L) {
    S <- Reduce(`+`, lapply(comps, function(z) make_full(z$ind, z$S)),
                init = matrix(0, p, p))
    return(list(S = list(S), lambda_mode = "scalar"))
  }

  if (q == n_sm) {
    S <- vector("list", n_sm)
    for (j in seq_len(n_sm)) S[[j]] <- matrix(0, p, p)
    for (z in comps) {
      S[[z$smooth]][z$ind, z$ind] <-
        S[[z$smooth]][z$ind, z$ind] + z$S
    }
    if (!is.null(setup$unpenalized) && setup$unpenalized > 0L) {
      ii <- seq_len(setup$unpenalized)
      S <- lapply(S, function(M) {
        M[ii, ] <- 0
        M[, ii] <- 0
        M
      })
    }
    return(list(S = S, lambda_mode = "smooth"))
  }

  if (q == n_pen) {
    S <- lapply(comps, function(z) make_full(z$ind, z$S))
    return(list(S = S, lambda_mode = "penalty"))
  }

  stop(
    "lambda_length must be 1, #smooths=", n_sm,
    ", or #penalties=", n_pen, "; got ", q, "."
  )
}
