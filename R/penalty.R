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

# Initial-fit penalty on the averaged score scale ------------------------------
#
# mgcv's smoothing parameters multiply the penalty on the *total* scale: gam()
# and bam() minimise the (unscaled) penalised deviance, or ||y - X beta||^2 +
# beta' S_lambda beta for the Gaussian, so the penalised Hessian of the
# log-likelihood is (X'WX + S_lambda) / phi and Vp = phi (X'WX + S_lambda)^-1,
# with W the IRLS weights computed from the unit variance function V(mu).
#
# The one-step equations in this package are on the *averaged* cluster scale,
# Wbar = N^-1 sum_i D_i' V_i^-1 D_i, where the working variance is
# V_i = c * V(mu) for a scalar c.  The pffr fit's own penalty on that scale is
# therefore S_lambda / (N * c).  This helper returns N * c:
#   * Gaussian: c = fit$sig2 (the working variance is the empirical residual
#     variance at each grid point, whose average is close to sig2);
#   * dispersion families: the dispersion multiplied into the working
#     variance (for example fit$sig2 for Gamma), and 1 if none was used, as for
#     binomial, Poisson and the quasi families with fixed nuisance;
#   * negative binomial and beta: 1, because theta/precision enter V(mu) in
#     the same way as in mgcv's extended families.
#' @keywords internal
#' @noRd
.fgee_initial_penalty_scale <- function(fit.initial, N, nuisance = NULL) {
  N <- as.numeric(N)[1L]
  if (!is.finite(N) || N < 1) stop("N must be a positive number of clusters.")
  key <- .fgee_family_key(fit.initial$family)
  c_scale <- 1
  if (key %in% c("gaussian", "normal")) {
    sig2 <- as.numeric(fit.initial$sig2)[1L]
    if (length(sig2) && is.finite(sig2) && sig2 > 0) c_scale <- sig2
  } else if (!key %in% c("negbinomial", "beta")) {
    disp <- if (is.list(nuisance)) as.numeric(nuisance$dispersion)[1L] else NA_real_
    if (length(disp) && is.finite(disp) && disp > 0) c_scale <- disp
  }
  N * c_scale
}
