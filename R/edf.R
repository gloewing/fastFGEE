#' Effective degrees of freedom by term (A_list) for fGEE
#' @keywords internal
#' @noRd
fgee_effective_df <- function(fit,
                              A_list,
                              wi_list = NULL,
                              penalty_mat = NULL,
                              A_col_tol = 0) {

  if (is.null(penalty_mat)) penalty_mat <- fit$pen.mat
  if (!is.list(A_list) || length(A_list) < 1) stop("A_list must be a non-empty list.")
  P <- as.matrix(penalty_mat)
  p <- ncol(P)

  if (is.null(wi_list)) wi_list <- fit$wi
  if (is.list(wi_list) && length(wi_list) > 0L) {
    Wbar <- Reduce(`+`, lapply(wi_list, as.matrix)) / length(wi_list)
  } else if (!is.null(fit$Wbar)) {
    Wbar <- as.matrix(fit$Wbar)
  } else if (inherits(fit$working, "fgee_working_stats")) {
    Wbar <- fit$working$W_bar
  } else {
    stop("Need wi_list, fit$Wbar, or fit$working.")
  }
  Wbar <- 0.5 * (Wbar + t(Wbar))
  P <- 0.5 * (P + t(P))

  S_mat <- solve_pd(Wbar + P, Wbar)
  dS <- diag(S_mat)

  R <- length(A_list)
  edf_r <- numeric(R)

  for (r in seq_len(R)) {
    Ar <- as.matrix(A_list[[r]])
    if (ncol(Ar) != p) stop("A_list[[", r, "]] has wrong ncol().")

    active <- which(colSums(abs(Ar)) > A_col_tol)
    if (!length(active)) active <- seq_len(p)

    edf_val <- sum(dS[active])
    edf_val <- max(0, min(edf_val, length(active)))
    edf_r[r] <- edf_val
  }

  list(edf_by_term = edf_r, diagS = dS)
}
