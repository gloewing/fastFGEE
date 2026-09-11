# Compact variance and wild-bootstrap helpers --------------------------------
#
# These functions operate on fgee_working_stats objects.  They avoid retaining
# one p x p W_i matrix per cluster for the default one-step sandwich and wild
# cluster bootstrap paths.

#' Variance estimation from compact working statistics
#' @keywords internal
#' @noRd
fgee_var_from_stats <- function(working,
                                beta,
                                penalty_diag,
                                working0 = NULL,
                                beta0 = NULL,
                                B = 500,
                                var.type = c("sandwich", "fastboot", "boot"),
                                exact = FALSE,
                                boot.base = c("initial", "final"),
                                return.boot = FALSE,
                                seed = NULL,
                                verbose = FALSE,
                                block_size = NULL,
                                max_counts_elems = 2e7) {
  if (!inherits(working, "fgee_working_stats")) {
    stop("working must be an fgee_working_stats object.")
  }
  var.type <- match.arg(var.type)
  boot.base <- match.arg(boot.base)
  P <- as.matrix(penalty_diag)
  beta <- as.numeric(beta)
  p <- working$p
  n <- working$N

  if (!all(dim(P) == c(p, p))) stop("penalty_diag must be p x p.")
  if (length(beta) != p) stop("beta must have length p.")
  if (!is.null(seed)) set.seed(seed)

  choose_block_size <- function(B, n_rows) {
    if (!is.null(block_size)) {
      bs <- as.integer(block_size)
      if (!is.finite(bs) || bs < 1L) stop("block_size must be positive.")
      return(min(bs, B))
    }
    max(1L, min(B, as.integer(floor(max_counts_elems / max(1L, n_rows)))))
  }

  if (var.type == "sandwich") {
    if (n < 2L) stop("var.type='sandwich' requires at least two clusters.")

    # The centered score crossproduct can be recovered from aggregate
    # sufficient statistics, so ordinary sandwich inference does not require
    # retaining the p x N score matrix.  Scores are still retained by default
    # for the studentized wild-cluster bootstrap.
    meat <- working$dd_sum - n * tcrossprod(working$d_bar)
    meat <- 0.5 * (meat + t(meat))

    A <- if (isTRUE(exact)) working$W_sum + P else working$W_bar + P
    fac <- .fgee_factor_pd(A)
    left <- fac$solve(meat)
    covhat <- t(fac$solve(t(left)))
    covhat <- 0.5 * (covhat + t(covhat))

    if (isTRUE(exact)) {
      return(covhat * (n / (n - 1L)))
    }
    return(covhat / (n * (n - 1L)))
  }

  if (var.type == "boot") {
    # True resample-and-resolve bootstrap changes the cluster-specific bread.
    # Delegate to the legacy implementation, which correctly requires W_i.
    W <- .fgee_working_wlist(working, required = TRUE)
    D <- .fgee_working_dlist(working, exact = isTRUE(exact), required = TRUE)
    W0 <- D0 <- NULL
    if (!is.null(working0)) {
      W0 <- .fgee_working_wlist(working0, required = TRUE)
      D0 <- .fgee_working_dlist(working0, required = TRUE)
    }
    return(var.est(
      di = D,
      wi = W,
      beta2 = beta,
      wi0 = W0,
      di0 = D0,
      beta0 = beta0,
      penalty_diag = P,
      B = B,
      var.type = "boot",
      exact = exact,
      boot.base = boot.base,
      return.boot = return.boot,
      seed = seed,
      verbose = verbose,
      block_size = block_size,
      max_counts_elems = max_counts_elems
    ))
  }

  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 2L) {
    stop("fastboot requires B >= 2.")
  }
  B <- as.integer(B)

  base <- working
  beta_start <- beta
  if (!isTRUE(exact) && boot.base == "initial" && !is.null(working0)) {
    base <- working0
    if (is.null(beta0)) stop("beta0 is required when boot.base='initial'.")
    beta_start <- as.numeric(beta0)
  }

  D <- .fgee_working_dmat(base, exact = isTRUE(exact), required = TRUE)
  n0 <- base$N
  bs <- choose_block_size(B, n0)
  prob <- rep.int(1 / n0, n0)

  if (!isTRUE(exact)) {
    Ainv <- solve_pd(base$W_bar + P)
    pen_term <- as.numeric(P %*% beta_start) * n0
  } else {
    Ainv <- solve_pd(base$W_sum + P)
    pen_term <- NULL
  }

  draws <- if (isTRUE(return.boot)) matrix(NA_real_, B, p) else NULL
  sum_x <- numeric(p)
  sum_xx <- matrix(0, p, p)

  b0 <- 1L
  while (b0 <= B) {
    Bb <- min(bs, B - b0 + 1L)
    counts <- stats::rmultinom(Bb, size = n0, prob = prob)
    Dsum <- D %*% counts

    if (!isTRUE(exact)) {
      beta_block <- Ainv %*% (Dsum - pen_term) / n0 + beta_start
    } else {
      beta_block <- Ainv %*% Dsum
    }

    if (!is.null(draws)) draws[b0:(b0 + Bb - 1L), ] <- t(beta_block)
    sum_x <- sum_x + rowSums(beta_block)
    sum_xx <- sum_xx + beta_block %*% t(beta_block)
    b0 <- b0 + Bb
  }

  covhat <- (sum_xx - tcrossprod(sum_x) / B) / (B - 1L)
  if (isTRUE(return.boot)) list(cov = covhat, boot = draws) else covhat
}

#' Effective degrees of freedom from an aggregate bread matrix
#' @keywords internal
#' @noRd
fgee_effective_df_from_Wbar <- function(Wbar,
                                        penalty_mat,
                                        A_list,
                                        A_col_tol = 0) {
  Wbar <- as.matrix(Wbar)
  P <- as.matrix(penalty_mat)
  if (!all(dim(Wbar) == dim(P))) stop("Wbar and penalty_mat must have equal dimensions.")
  if (!is.list(A_list) || !length(A_list)) stop("A_list must be a non-empty list.")

  Wbar <- 0.5 * (Wbar + t(Wbar))
  P <- 0.5 * (P + t(P))
  S <- solve_pd(Wbar + P, Wbar)
  dS <- diag(S)
  p <- ncol(P)

  edf <- vapply(seq_along(A_list), function(r) {
    A <- as.matrix(A_list[[r]])
    if (ncol(A) != p) stop("A_list[[", r, "]] has the wrong number of columns.")
    active <- which(colSums(abs(A)) > A_col_tol)
    if (!length(active)) active <- seq_len(p)
    max(0, min(sum(dS[active]), length(active)))
  }, numeric(1))

  list(edf_by_term = edf, diagS = dS)
}

#' Studentized wild cluster bootstrap from compact statistics
#' @keywords internal
#' @noRd
wild_studentized_boot_crit_compact <- function(beta_hat,
                                                Sigma_hat = NULL,
                                                D,
                                                Wbar,
                                                penalty_diag,
                                                A_list,
                                                B = 5000,
                                                seed = NULL,
                                                alpha = 0.05,
                                                studentize = c("fixed", "replicate"),
                                                se_floor_rel = 1e-6,
                                                progress_every = 1000,
                                                return_T = FALSE,
                                                return_bias = FALSE,
                                                beta_center = NULL,
                                                t_adjust = c("both", "joint", "none", "pointwise"),
                                                df_floor = 1,
                                                edf_r = NULL,
                                                A_col_tol = 0,
                                                return_ci = TRUE) {
  studentize <- match.arg(studentize)
  t_adjust <- match.arg(t_adjust)
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0,1).")
  }

  D <- as.matrix(D)
  Wbar <- as.matrix(Wbar)
  P <- as.matrix(penalty_diag)
  beta_hat <- as.numeric(beta_hat)
  p <- length(beta_hat)
  n <- ncol(D)

  if (nrow(D) != p) stop("D must be p x N.")
  if (!all(dim(Wbar) == c(p, p))) stop("Wbar must be p x p.")
  if (!all(dim(P) == c(p, p))) stop("penalty_diag must be p x p.")
  if (!is.list(A_list) || !length(A_list)) stop("A_list must be non-empty.")
  if (B < 2L) stop("B must be at least 2.")

  if (is.null(beta_center)) beta_center <- beta_hat
  beta_center <- as.numeric(beta_center)
  if (length(beta_center) != p) stop("beta_center must have length p.")

  P <- 0.5 * (P + t(P))
  Wbar <- 0.5 * (Wbar + t(Wbar))
  Hinv <- solve_pd(Wbar + P)
  penvec <- as.numeric(P %*% beta_hat) * n

  if (is.null(Sigma_hat)) {
    Sigma_hat <- Hinv %*% tcrossprod(D) %*% Hinv / (n * (n - 1L))
  } else {
    Sigma_hat <- as.matrix(Sigma_hat)
    if (!all(dim(Sigma_hat) == c(p, p))) stop("Sigma_hat must be p x p.")
  }

  R <- length(A_list)
  nr <- vapply(A_list, nrow, integer(1))
  if (length(unique(nr)) != 1L) stop("All A_list matrices must have the same nrow().")
  m <- nr[1L]

  fhat_pivot <- sehat <- vector("list", R)
  for (r in seq_len(R)) {
    A <- as.matrix(A_list[[r]])
    if (ncol(A) != p) stop("A_list[[", r, "]] must have p columns.")
    fhat_pivot[[r]] <- as.numeric(A %*% beta_hat)
    se <- sqrt(pmax(rowSums((A %*% Sigma_hat) * A), 0))
    sehat[[r]] <- pmax(se, se_floor_rel * se + .Machine$double.eps)
  }

  T_list <- lapply(seq_len(R), function(i) matrix(NA_real_, B, m))
  M_mat <- matrix(NA_real_, B, R)
  f_sum <- if (isTRUE(return_bias)) lapply(seq_len(R), function(i) numeric(m)) else NULL

  for (b in seq_len(B)) {
    xi <- sample(c(-1, 1), n, replace = TRUE)
    beta_b <- beta_hat + as.numeric(Hinv %*% (D %*% xi - penvec)) / n

    for (r in seq_len(R)) {
      A <- as.matrix(A_list[[r]])
      f_b <- as.numeric(A %*% beta_b)
      se_b <- sehat[[r]]
      if (studentize == "replicate") {
        # The current package's replicate branch also keeps Sigma fixed.
        se_b <- sqrt(pmax(rowSums((A %*% Sigma_hat) * A), 0))
        se_b <- pmax(se_b, se_floor_rel * sehat[[r]])
      }
      if (isTRUE(return_bias)) f_sum[[r]] <- f_sum[[r]] + f_b
      tr <- (f_b - fhat_pivot[[r]]) / se_b
      T_list[[r]][b, ] <- tr
      M_mat[b, r] <- max(abs(tr), na.rm = TRUE)
    }

    if (progress_every > 0L && b %% progress_every == 0L) {
      message("wild bootstrap iter: ", b, " / ", B)
    }
  }

  crit <- list(
    joint = vector("list", R),
    pt_lo = vector("list", R),
    pt_hi = vector("list", R),
    pointwise_crit = vector("list", R),
    info = list(studentize = studentize, B = B, alpha = alpha)
  )
  for (r in seq_len(R)) {
    tr <- T_list[[r]]
    crit$pt_lo[[r]] <- apply(tr, 2L, stats::quantile, probs = alpha / 2, na.rm = TRUE)
    crit$pt_hi[[r]] <- apply(tr, 2L, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    crit$joint[[r]] <- as.numeric(stats::quantile(M_mat[, r], 1 - alpha, na.rm = TRUE))
    crit$pointwise_crit[[r]] <- as.numeric(stats::quantile(abs(as.vector(tr)), 1 - alpha, na.rm = TRUE))
  }

  if (is.null(edf_r)) {
    edf_r <- fgee_effective_df_from_Wbar(
      Wbar = Wbar,
      penalty_mat = P,
      A_list = A_list,
      A_col_tol = A_col_tol
    )$edf_by_term
  }
  edf_r <- as.numeric(edf_r)
  if (length(edf_r) != R) stop("edf_r must have one value per A_list element.")
  df_r <- pmax(df_floor, n - edf_r)
  adj_r <- stats::qt(1 - alpha / 2, df_r) / stats::qnorm(1 - alpha / 2)

  crit$info$edf_r <- edf_r
  crit$info$df_r <- df_r
  crit$info$adj_r <- adj_r
  crit$info$t_adjust <- t_adjust
  crit$info$beta_centered_at <- beta_center
  crit$joint_adj <- crit$pointwise_crit_adj <- vector("list", R)

  for (r in seq_len(R)) {
    jc <- crit$joint[[r]]
    pc <- crit$pointwise_crit[[r]]
    if (t_adjust %in% c("joint", "both")) jc <- jc * adj_r[r]
    if (t_adjust %in% c("pointwise", "both")) pc <- pc * adj_r[r]
    crit$joint_adj[[r]] <- as.numeric(jc)
    crit$pointwise_crit_adj[[r]] <- as.numeric(pc)
  }

  if (isTRUE(return_ci)) {
    nm <- names(A_list)
    if (is.null(nm)) nm <- paste0("term", seq_len(R))
    fit_center <- ci_pointwise <- ci_joint <- df_list <- vector("list", R)

    for (r in seq_len(R)) {
      A <- as.matrix(A_list[[r]])
      fit_r <- as.numeric(A %*% beta_center)
      se_r <- sehat[[r]]
      cp <- crit$pointwise_crit_adj[[r]]
      cj <- crit$joint_adj[[r]]
      lo_p <- fit_r - cp * se_r
      hi_p <- fit_r + cp * se_r
      lo_j <- fit_r - cj * se_r
      hi_j <- fit_r + cj * se_r
      fit_center[[r]] <- fit_r
      ci_pointwise[[r]] <- cbind(lower = lo_p, upper = hi_p)
      ci_joint[[r]] <- cbind(lower = lo_j, upper = hi_j)
      df_list[[r]] <- data.frame(
        term = nm[r],
        s = seq_len(length(fit_r)),
        fit = fit_r,
        se = se_r,
        lower_pt = lo_p,
        upper_pt = hi_p,
        lower_joint = lo_j,
        upper_joint = hi_j,
        edf = edf_r[r],
        df = df_r[r],
        adj = adj_r[r],
        pointwise_crit = cp,
        joint_crit = cj
      )
    }
    names(fit_center) <- names(ci_pointwise) <- names(ci_joint) <- names(df_list) <- nm
    crit$ci <- list(
      fit = fit_center,
      se = sehat,
      ci_pointwise = ci_pointwise,
      ci_joint = ci_joint,
      df = df_list,
      info = list(
        alpha = alpha,
        t_adjust = t_adjust,
        centered_at = "beta_center",
        pivot_center = "beta_hat"
      ),
      fit_pivot = fhat_pivot
    )
  }

  if (isTRUE(return_T)) {
    attr(crit, "T_list") <- T_list
    attr(crit, "M_mat") <- M_mat
  }
  if (isTRUE(return_bias)) crit$bias_sum <- f_sum
  crit
}

#' Compact joint-CI dispatcher
#' @keywords internal
#' @noRd
compute_joint_ci_compact <- function(joint.CI = FALSE,
                                     glmfit,
                                     working0,
                                     working2,
                                     beta0,
                                     beta2,
                                     penalty_diag,
                                     bs = NULL,
                                     index = "yindex.vec",
                                     alpha = 0.05,
                                     wild_B = 5000,
                                     wild_seed = 123,
                                     wild_progress_every = 1000,
                                     t_adjust = c("both", "joint", "none", "pointwise"),
                                     df_floor = 1,
                                     se_floor_rel = 1e-6,
                                     A_col_tol = 0,
                                     exact = FALSE) {
  t_adjust <- match.arg(t_adjust)
  if (identical(joint.CI, FALSE)) return(list(qn = NULL, glmfit = glmfit))

  if (!identical(joint.CI, "wild")) {
    # Preserve every legacy CI type when full cluster matrices were requested.
    return(compute_joint_ci(
      joint.CI = joint.CI,
      glmfit = glmfit,
      di0 = .fgee_working_dlist(working0, required = TRUE),
      wi0 = .fgee_working_wlist(working0, required = TRUE),
      di2 = .fgee_working_dlist(working2, required = TRUE),
      wi2 = .fgee_working_wlist(working2, required = TRUE),
      beta0 = beta0,
      beta2 = beta2,
      penalty_diag = penalty_diag,
      bs = bs,
      index = index,
      alpha = alpha,
      t_adjust = t_adjust,
      df_floor = df_floor,
      se_floor_rel = se_floor_rel,
      A_col_tol = A_col_tol,
      exact = exact
    ))
  }

  message("CIs: studentized wild cluster bootstrap (compact pointwise + joint)")
  dd_ci <- make_ci_newdata(glmfit, fn_domain = index, grid.size = NULL)
  dd_ci <- droplevels(dd_ci)
  by_vars <- unique(stats::na.omit(sapply(glmfit$smooth, function(s) s$by)))
  by_vars <- intersect(by_vars, names(dd_ci))
  for (bv in by_vars) if (is.numeric(dd_ci[[bv]])) dd_ci[[bv]] <- 1
  A_list <- build_A_functional_by(glmfit, newdata = dd_ci, fn_domain = index)

  edf <- fgee_effective_df_from_Wbar(
    Wbar = working2$W_bar,
    penalty_mat = penalty_diag,
    A_list = A_list,
    A_col_tol = A_col_tol
  )

  crit <- wild_studentized_boot_crit_compact(
    beta_hat = as.numeric(beta0),
    Sigma_hat = glmfit$Vp,
    D = .fgee_working_dmat(working0, required = TRUE),
    Wbar = working0$W_bar,
    penalty_diag = penalty_diag,
    A_list = A_list,
    B = wild_B,
    seed = wild_seed,
    alpha = alpha,
    return_T = FALSE,
    beta_center = as.numeric(beta2),
    t_adjust = if (isTRUE(exact)) "none" else t_adjust,
    edf_r = edf$edf_by_term,
    studentize = "fixed",
    progress_every = wild_progress_every,
    df_floor = df_floor,
    se_floor_rel = se_floor_rel,
    A_col_tol = A_col_tol
  )

  glmfit$crit <- crit
  glmfit$ci_newdata <- dd_ci
  glmfit$ci <- crit$ci
  glmfit$ci_fn_domain <- index
  glmfit$ci_names <- names(A_list)
  list(qn = NULL, glmfit = glmfit)
}
