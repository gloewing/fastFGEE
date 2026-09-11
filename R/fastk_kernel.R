# Optional compiled kernel for the fastK objective ---------------------------
#
# The per-fold work in `fgee_fastk_score_grad()` is, for each fold, a matrix
# slice, a matrix-vector product for eta, a vectorised loss and derivative, and
# a crossproduct for X'g.  That allocates three or four vectors of length M_k per
# fold per evaluation.  The compiled kernel does the same arithmetic in a single
# pass with no intermediate allocation.
#
# The kernel is an accelerator, not a change of method.  It reproduces the R
# objective to about 1e-16 in both the score and the gradient, so enabling or
# disabling it cannot move a fitted model.  Every entry point degrades to the R
# implementation when the compiled code is unavailable, which keeps the package
# usable if the shared object fails to load for any reason.

#' Is the compiled fastK kernel available?
#'
#' @return `TRUE` when the compiled loss/gradient kernel can be used.
#' @keywords internal
#' @noRd
fgee_fastk_kernel_ok <- function() {
  isTRUE(getOption("fastFGEE.kernel", TRUE)) &&
    is.function(tryCatch(get(".fastk_fold_kernel", envir = asNamespace("fastFGEE")),
                         error = function(e) NULL))
}

#' Map a family key to the kernel's integer family code
#' @keywords internal
#' @noRd
.fgee_fastk_family_code <- function(key) {
  switch(key,
    gaussian = 0L, normal = 0L,
    binomial = 1L, quasibinomial = 1L,
    poisson = 2L, quasipoisson = 2L,
    gamma = 3L, quasigamma = 3L,
    NULL)
}

#' Exact fastK score and gradient via the compiled kernel
#'
#' Identical in arithmetic to [fgee_fastk_score_grad()]; only the per-fold loss
#' and `X'g` accumulation is delegated.  Requires `memory = "balanced"`, which
#' supplies the single contiguous design matrix the kernel reads.
#'
#' @keywords internal
#' @noRd
.fgee_fastk_score_grad_kernel <- function(prep, lambda, need_gradient = TRUE,
                                          return_fold = FALSE) {
  lambda <- as.numeric(lambda)
  if (length(lambda) != prep$q || any(!is.finite(lambda)) || any(lambda <= 0)) {
    return(list(score = Inf, gradient = rep(NA_real_, prep$q)))
  }
  fam <- .fgee_fastk_family_code(prep$family_key)
  if (is.null(prep$X_eval) || is.null(fam) || isTRUE(prep$gaussian)) {
    return(fgee_fastk_score_grad(prep, lambda, need_gradient = need_gradient,
                                 return_fold = return_fold))
  }

  P <- penalty_from_setup(prep$ps, lambda)
  if (isTRUE(prep$exact)) {
    fac <- .fgee_factor_pd(prep$Wbar + P / prep$ncl)
    Bhat <- fac$solve(prep$D / prep$ncl)
  } else {
    fac <- .fgee_factor_pd(prep$Wbar + P)
    Bhat <- prep$beta_mat +
      fac$solve(prep$D - (P %*% prep$beta_mat) * prep$ncl) / prep$ncl
  }

  kk <- .fastk_fold_kernel(prep$X_eval, prep$Y_eval, prep$W_eval,
                           prep$fold_starts, prep$fold_ends, Bhat,
                           fam, prep$clip_prob, need_gradient)

  score <- sum(kk$loss)
  fold_weight <- vapply(prep$folds$sets, length, integer(1)) / prep$ncl
  fold_mean <- kk$loss / fold_weight

  if (!need_gradient) {
    ans <- list(score = score, beta_folds = Bhat)
    if (return_fold) ans$fold_score <- fold_mean
    return(ans)
  }

  grad <- numeric(prep$q)
  l10 <- log(10)
  for (j in seq_len(prep$q)) {
    dP <- (l10 * lambda[j]) * prep$components[[j]]
    if (isTRUE(prep$exact)) dP <- dP / prep$ncl
    grad[j] <- sum(kk$G * (-fac$solve(dP %*% Bhat)))
  }
  ans <- list(score = score, gradient = grad, beta_folds = Bhat, Gbeta = kk$G)
  if (return_fold) ans$fold_score <- fold_mean
  ans
}
