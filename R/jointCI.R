# Helper function for computing joint confidence intervals
#' @keywords internal
#' @noRd
compute_joint_ci <- function(
    joint.CI = FALSE,
    glmfit,
    di0,
    wi0,
    di2,
    wi2,
    beta0,
    beta2,
    penalty_diag,
    bs = NULL,
    index = "yindex.vec",
    alpha = 0.05,
    parametric_nsim = 10000,
    np_nsim = 10000,
    np_grid_size = 500,
    wild_B = 5000,
    wild_seed = 123,
    wild_progress_every = 1000,
    t_adjust = c("both", "joint", "none", "pointwise"),
    df_floor = 1,
    se_floor_rel = 1e-6,
    A_col_tol = 0,
    exact = FALSE
) {

  t_adjust <- match.arg(t_adjust)

  qn <- NULL

  if (identical(joint.CI, FALSE)) {
    return(list(qn = NULL, glmfit = glmfit))
  }

  # ----------------------------
  # NON-WILD joint criticals (qn)
  # ----------------------------
  if (joint.CI == "parametric") {

    message("joint CIs: parametric bootstrap (critical values)")
    qn <- joint.qn(glmfit = glmfit, NN.sim = parametric_nsim, alpha = alpha)

    # Attach wild-shaped CI object using Wald + qn
    out <- .ci_build_wald(
      glmfit = glmfit,
      wi_list = wi2,
      penalty_diag = penalty_diag,
      beta_center = as.numeric(beta2),
      fn_domain = index,
      alpha = alpha,
      qn_joint = qn,
      qn_label = "joint.qn(parametric)",
      t_adjust = t_adjust,
      df_floor = df_floor,
      se_floor_rel = se_floor_rel,
      A_col_tol = A_col_tol
    )

    glmfit$crit <- out$crit
    glmfit$ci_newdata <- out$dd_ci
    glmfit$ci <- out$crit$ci
    glmfit$ci_fn_domain <- index
    glmfit$ci_names <- names(out$A_list)

    # IMPORTANT: Return the updated glmfit with qn
    return(list(qn = qn, glmfit = glmfit))

  } else if (joint.CI == "np") {

    message("joint CIs: non-parametric bootstrap (basis / non-studentized critical values)")
    qn <- joint.basis.np(
      glmfit = glmfit,
      di = di2,
      wi = wi2,
      beta = beta2,
      penalty_diag = penalty_diag,
      bs = bs,
      var.mat = glmfit$Vp,
      exact = FALSE,
      grid.size = np_grid_size,
      NN.sim = np_nsim,
      alpha = alpha
    )

    out <- .ci_build_wald(
      glmfit = glmfit,
      wi_list = wi2,
      penalty_diag = penalty_diag,
      beta_center = as.numeric(beta2),
      fn_domain = index,
      alpha = alpha,
      qn_joint = qn,
      qn_label = "joint.basis.np(np)",
      t_adjust = t_adjust,
      df_floor = df_floor,
      se_floor_rel = se_floor_rel,
      A_col_tol = A_col_tol
    )

    glmfit$crit <- out$crit
    glmfit$ci_newdata <- out$dd_ci
    glmfit$ci <- out$crit$ci
    glmfit$ci_fn_domain <- index
    glmfit$ci_names <- names(out$A_list)

    # IMPORTANT: Return the updated glmfit with qn
    return(list(qn = qn, glmfit = glmfit))

  } else if (joint.CI == "np_curve") {

    # Optional: uses your existing joint.np() (grid-based curve criticals)
    message("joint CIs: non-parametric bootstrap (curve / grid critical values)")
    qn <- joint.np(
      glmfit = glmfit,
      di = di2,
      wi = wi2,
      beta = beta2,
      penalty_diag = penalty_diag,
      bs = bs,
      var.mat = glmfit$Vp,
      grid.size = np_grid_size,
      NN.sim = np_nsim,
      alpha = alpha
    )

    out <- .ci_build_wald(
      glmfit = glmfit,
      wi_list = wi2,
      penalty_diag = penalty_diag,
      beta_center = as.numeric(beta2),
      fn_domain = index,
      alpha = alpha,
      qn_joint = qn,
      qn_label = "joint.np(np_curve)",
      t_adjust = t_adjust,
      df_floor = df_floor,
      se_floor_rel = se_floor_rel,
      A_col_tol = A_col_tol
    )

    glmfit$crit <- out$crit
    glmfit$ci_newdata <- out$dd_ci
    glmfit$ci <- out$crit$ci
    glmfit$ci_fn_domain <- index
    glmfit$ci_names <- names(out$A_list)

    # IMPORTANT: Return the updated glmfit with qn
    return(list(qn = qn, glmfit = glmfit))

  } else if (joint.CI == "wild") {

    # ----------------------------
    # WILD BRANCH (UNCHANGED)
    # ----------------------------
    message("CIs: studentized wild cluster bootstrap (pointwise + joint)")

    dd_ci <- make_ci_newdata(glmfit, fn_domain = index, grid.size = NULL)
    dd_ci <- droplevels(dd_ci)

    by_vars <- unique(na.omit(sapply(glmfit$smooth, function(s) s$by)))
    by_vars <- intersect(by_vars, names(dd_ci))
    for (bv in by_vars) {
      if (is.numeric(dd_ci[[bv]])) dd_ci[[bv]] <- 1
    }

    A_list <- build_A_functional_by(glmfit, newdata = dd_ci, fn_domain = index)

    edf_r <- fgee_effective_df(
      glmfit,
      A_list = A_list,
      wi_list = wi2,
      penalty_mat = penalty_diag
    )

    crit_wild <- wild_studentized_boot_crit(
      beta_hat = as.numeric(beta0), # if (exact) as.numeric(beta2) else as.numeric(beta0),
      Sigma_hat = glmfit$Vp,
      di = di0, #if (exact) di2 else di0, # di0 is correct even for exact (calculated as in for 1-step)
      wi = wi0, #if (exact) wi2 else wi0,
      penalty_diag = penalty_diag,
      A_list = A_list,
      B = wild_B,
      seed = wild_seed,
      return_T = FALSE,
      beta_center = as.numeric(beta2),
      t_adjust = if (exact) "none" else "both", #
      edf_r = edf_r$edf_by_term,
      studentize = "fixed",
      progress_every = wild_progress_every
    )

    nr <- vapply(A_list, nrow, integer(1))
    if (length(unique(nr)) != 1) stop("A_list matrices have inconsistent nrow().")

    glmfit$crit <- crit_wild
    glmfit$ci_newdata <- dd_ci
    glmfit$ci <- crit_wild$ci
    glmfit$ci_fn_domain <- index
    glmfit$ci_names <- names(A_list)

    qn <- NULL

    # Return for wild bootstrap
    return(list(qn = qn, glmfit = glmfit))

  } else {
    stop("Unsupported joint.CI='", joint.CI, "'. Supported: FALSE, 'parametric', 'np', 'np_curve', 'wild'.")
  }
}
#' @keywords internal
#' @noRd
joint.qn <- function(glmfit, NN.sim = 5000, alpha = 0.05) {
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }

  p.vec <- sapply(glmfit$smooth, function(xx) ncol(xx$S[[1]]))
  p.vec[1] <- p.vec[1] + 1  # your original intercept adjustment
  p.vec <- c(0, cumsum(p.vec))

  qn <- rep(0, length = length(glmfit$smooth))

  for (i in seq_along(qn)) {
    cov.idx <- (p.vec[i] + 1):(p.vec[i + 1])
    Sigma <- glmfit$Vp[cov.idx, cov.idx, drop = FALSE]

    sqrt_Sigma <- sqrt(diag(Sigma))
    S_scl <- Matrix::Diagonal(x = 1 / sqrt_Sigma)
    Sigma_scl <- as.matrix(S_scl %*% Sigma %*% S_scl)

    zero_vec <- rep(0, nrow(Sigma_scl))
    x_sample <- abs(MASS::mvrnorm(NN.sim, zero_vec, Sigma_scl))
    un <- Rfast::rowMaxs(x_sample, value = TRUE)
    qn[i] <- stats::quantile(un, probs = 1 - alpha)
  }

  qn
}

#' @keywords internal
#' @noRd
joint.basis.np <- function(glmfit,
                           di,
                           wi,
                           beta,
                           penalty_diag,
                           bs,
                           var.mat,
                           exact = FALSE,
                           grid.size = 500,
                           NN.sim = 5000,
                           alpha = 0.05) {

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }

  boot.samp <- suppressMessages(
    var.est(di = di, wi = wi, beta = beta,
            penalty_diag = penalty_diag, B = NN.sim,
            var.type = "fastboot", exact = exact,
            return.boot = TRUE)
  )$boot

  pp <- length(glmfit$smooth)
  idx.cs <- sapply(seq_len(pp), function(p) ncol(glmfit$smooth[[p]]$D))
  idx.cs <- c(0, cumsum(idx.cs))

  qn <- rep(NA_real_, pp)
  var.vec <- diag(var.mat)

  for (k in seq_len(pp)) {
    b.idx <- seq(idx.cs[k] + 1, idx.cs[k + 1])
    b.mat <- boot.samp[, b.idx, drop = FALSE]

    z.beta <- abs(beta[b.idx] - t(b.mat)) / sqrt(var.vec[b.idx])
    un <- Rfast::colMaxs(z.beta, value = TRUE)
    qn[k] <- stats::quantile(un, probs = 1 - alpha)
  }

  qn
}

#' @keywords internal
#' @noRd
joint.np <- function(glmfit,
                     di,
                     wi,
                     beta,
                     penalty_diag,
                     bs,
                     var.mat,
                     grid.size = 500,
                     NN.sim = 5000,
                     alpha = 0.05) {

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }

  boot.samp <- suppressMessages(
    var.est(di = di, wi = wi, beta = beta,
            penalty_diag = penalty_diag, B = NN.sim,
            var.type = "fastboot", exact = FALSE,
            return.boot = TRUE)
  )
  var.mat <- boot.samp$cov
  boot.samp <- boot.samp$boot

  # NOTE: your original code uses glmfit$pffr$yind; keep as-is for now.
  yind <- glmfit$pffr$yind
  if (length(yind) > grid.size) grid.size <- length(yind)

  x <- seq(min(yind), max(yind), length = grid.size)

  pp <- length(glmfit$smooth)
  idx <- sapply(seq_len(pp), function(p) ncol(glmfit$smooth[[p]]$D))
  idx.cs <- c(0, cumsum(idx))
  grid.idx <- c(0, cumsum(rep(grid.size, pp)))

  qn <- rep(NA_real_, pp)

  sm.X.ls <- lapply(seq_len(pp), function(kk)
    mgcv::smoothCon(s(x, bs = bs, k = idx[kk]), data = data.frame(x))[[1]]$X
  )

  XX <- Matrix::bdiag(sm.X.ls)
  var.vec <- rowSums(XX * t(var.mat %*% t(XX)))
  rm(XX)

  for (k in seq_len(pp)) {
    b.idx <- seq(idx.cs[k] + 1, idx.cs[k + 1])
    b.mat <- boot.samp[, b.idx, drop = FALSE]
    sm.X <- sm.X.ls[[k]]

    f.boot.hat <- tcrossprod(sm.X, b.mat)

    grid.idx.k <- seq(grid.idx[k] + 1, grid.idx[k + 1])
    f.hat <- as.numeric(sm.X %*% beta[b.idx])

    z.beta <- abs(f.hat - f.boot.hat) / sqrt(var.vec[grid.idx.k])
    un <- Rfast::colMaxs(z.beta, value = TRUE)
    qn[k] <- stats::quantile(un, probs = 1 - alpha)
  }

  qn
}




#' @keywords internal
#' @noRd
wild_studentized_boot_crit <- function(beta_hat, Sigma_hat = NULL,
                                       di, wi, penalty_diag,
                                       A_list,
                                       B = 5000, seed = NULL,
                                       alpha = 0.05,
                                       studentize = c("fixed", "replicate"),
                                       se_floor_rel = 1e-6,
                                       progress_every = 1000,
                                       return_T = FALSE,
                                       return_bias = FALSE,
                                       beta_center = NULL, # <-- ONLY affects CI centering
                                       # ---------------- NEW ----------------
                                       t_adjust = c("both", "joint", "none", "pointwise"),
                                       df_floor = 1,
                                       edf_r = NULL,
                                       A_col_tol = 0,
                                       return_ci = TRUE) {

  studentize <- match.arg(studentize)
  t_adjust   <- match.arg(t_adjust)

  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0,1).")
  }

  stopifnot(is.list(di), is.list(wi), is.list(A_list))
  n <- length(di)
  if (length(wi) != n) stop("di and wi must have the same length (#clusters).")
  if (B < 2) stop("Need B >= 2.")

  beta_hat <- as.numeric(beta_hat)
  p <- length(beta_hat)

  # CI-centering beta (does NOT affect bootstrap pivots/criticals)
  if (is.null(beta_center)) {
    beta_center <- beta_hat
  } else {
    beta_center <- as.numeric(beta_center)
    if (length(beta_center) != p) stop("beta_center must have length = length(beta_hat).")
  }

  # penalty matrix (tuned, not mgcv), symmetrize for stability
  P <- as.matrix(penalty_diag)
  if (!all(dim(P) == c(p, p))) stop("penalty_diag must be p x p.")
  P <- 0.5 * (P + t(P))

  # stack di into p x n matrix for fast signed sums
  U_mat <- do.call(cbind, lapply(di, function(u) {
    u <- as.numeric(u)
    if (length(u) != p) stop("Each di[[i]] must have length p.")
    u
  }))  # p x n

  # bread: Wbar + P
  Wbar <- Reduce("+", lapply(wi, function(Wi) {
    Wi <- as.matrix(Wi)
    if (!all(dim(Wi) == c(p, p))) stop("Each wi[[i]] must be p x p.")
    Wi
  })) / n
  Wbar <- 0.5 * (Wbar + t(Wbar))

  H <- Wbar + P

  # invert bread
  Winv <- solve(H)

  # penalty vector for one-step (same scaling as your old code)
  penvec <- as.numeric(P %*% beta_hat * n)

  # Sigma_hat fallback if NULL: 1-step sandwich-like (same as before)
  if (is.null(Sigma_hat)) {
    U_orig <- U_mat %*% t(U_mat)  # sum_i di_i di_i^T
    Sigma_hat <- (1 / (n * (n - 1))) * (Winv %*% U_orig %*% Winv)
  } else {
    Sigma_hat <- as.matrix(Sigma_hat)
    if (!all(dim(Sigma_hat) == c(p, p))) stop("Sigma_hat must be p x p.")
  }

  # Precompute pivot-centering fhat (beta_hat) + SEs (Sigma_hat)
  # NOTE: these DO affect pivots and therefore critical values -> keep at beta_hat
  R <- length(A_list)
  if (R < 1) stop("A_list must be non-empty.")
  m_vec <- vapply(A_list, nrow, integer(1))
  if (length(unique(m_vec)) != 1) stop("All A_list matrices must have same nrow().")
  m <- m_vec[1]

  fhat_pivot  <- vector("list", R)  # A_r %*% beta_hat (used inside pivots)
  sehat       <- vector("list", R)  # sqrt(diag(A Sigma A^T)) (used for widths)
  for (r in seq_len(R)) {
    Ar <- as.matrix(A_list[[r]])
    if (ncol(Ar) != p) stop("A_list[[", r, "]] must have ncol = p.")

    fhat_pivot[[r]] <- as.numeric(Ar %*% beta_hat)

    se <- sqrt(pmax(diag(Ar %*% Sigma_hat %*% t(Ar)), 0))
    se <- pmax(se, se_floor_rel * se + .Machine$double.eps)
    sehat[[r]] <- se
  }

  # storage
  T_list <- lapply(seq_len(R), function(i) matrix(NA_real_, nrow = B, ncol = m))
  M_mat  <- matrix(NA_real_, nrow = B, ncol = R)
  f_sum  <- if (isTRUE(return_bias)) lapply(seq_len(R), function(i) numeric(m)) else NULL

  # main loop (bootstrap run unchanged)
  for (b in seq_len(B)) {
    xi <- sample(c(-1, 1), n, replace = TRUE)   # Rademacher
    Dsum_b <- as.numeric(U_mat %*% xi)          # signed sum (p-vector)

    # one-step update (fixed bread)
    inf_b  <- Winv %*% (Dsum_b - penvec)
    beta_b <- beta_hat + as.numeric(inf_b) / n

    # replicate SE option (kept consistent with your prior "replicate" logic)
    if (studentize == "replicate") {
      Sig_b <- Sigma_hat
    } else {
      Sig_b <- NULL
    }

    for (r in seq_len(R)) {
      Ar <- as.matrix(A_list[[r]])
      f_b <- as.numeric(Ar %*% beta_b)

      if (studentize == "replicate") {
        se_b <- sqrt(pmax(diag(Ar %*% Sig_b %*% t(Ar)), 0))
        se_b <- pmax(se_b, se_floor_rel * sehat[[r]])
      } else {
        se_b <- sehat[[r]]
      }

      if (isTRUE(return_bias)) f_sum[[r]] <- f_sum[[r]] + f_b

      # pivots centered at fhat_pivot (beta_hat) — unchanged
      Tr <- (f_b - fhat_pivot[[r]]) / se_b
      T_list[[r]][b, ] <- Tr
      M_mat[b, r] <- max(abs(Tr), na.rm = TRUE)
    }

    if (progress_every > 0 && (b %% progress_every == 0)) {
      message("wild bootstrap iter: ", b, " / ", B)
    }
  }

  # ------------------------------------------------------------
  # critical values (unchanged) + scalar pointwise crit (Fix 1)
  # ------------------------------------------------------------
  crit <- list(
    joint = vector("list", R),
    pt_lo = vector("list", R),
    pt_hi = vector("list", R),
    pointwise_crit = vector("list", R),
    info = list(studentize = studentize, B = B, alpha = alpha)
  )

  for (r in seq_len(R)) {
    Tr <- T_list[[r]]
    crit$pt_lo[[r]] <- apply(Tr, 2, stats::quantile, probs = alpha/2, na.rm = TRUE)
    crit$pt_hi[[r]] <- apply(Tr, 2, stats::quantile, probs = 1 - alpha/2, na.rm = TRUE)
    crit$joint[[r]] <- as.numeric(stats::quantile(M_mat[, r], probs = 1 - alpha, na.rm = TRUE))

    # EXACTLY your Fix 1 scalar pointwise critical value:
    crit$pointwise_crit[[r]] <- as.numeric(stats::quantile(abs(as.vector(Tr)),
                                                           probs = 1 - alpha, na.rm = TRUE))
  }

  # ------------------------------------------------------------
  # df-based adjustment using tuned penalty matrix P
  # ------------------------------------------------------------
  if (is.null(edf_r)) {
    S_mat <- Winv %*% Wbar
    dS <- diag(S_mat)

    edf_r <- numeric(R)
    for (r in seq_len(R)) {
      Ar <- as.matrix(A_list[[r]])
      active <- which(colSums(abs(Ar)) > A_col_tol)
      if (!length(active)) active <- seq_len(p)
      edf_val <- sum(dS[active])
      edf_val <- max(0, min(edf_val, length(active)))
      edf_r[r] <- edf_val
    }
  } else {
    edf_r <- as.numeric(edf_r)
    if (length(edf_r) != R) stop("edf_r must have length = length(A_list).")
  }

  df_r <- pmax(df_floor, n - edf_r)
  adj_r <- stats::qt(1 - alpha/2, df = df_r) / stats::qnorm(1 - alpha/2)

  crit$info$edf_r <- edf_r
  crit$info$df_r  <- df_r
  crit$info$adj_r <- adj_r
  crit$info$t_adjust <- t_adjust
  crit$info$beta_centered_at <- beta_center

  # adjusted critical values
  crit$joint_adj <- vector("list", R)
  crit$pointwise_crit_adj <- vector("list", R)
  for (r in seq_len(R)) {
    jcrit <- as.numeric(crit$joint[[r]])
    pcrit <- as.numeric(crit$pointwise_crit[[r]])

    if (t_adjust %in% c("joint", "both")) {
      jcrit <- jcrit * adj_r[r]
    }
    if (t_adjust %in% c("pointwise", "both")) {
      pcrit <- pcrit * adj_r[r]
    }

    crit$joint_adj[[r]] <- jcrit
    crit$pointwise_crit_adj[[r]] <- pcrit
  }

  # ------------------------------------------------------------
  # Construct CIs using:
  #   - widths from bootstrap (criticals + SEs)
  #   - centers from beta_center (THIS is the only change you requested)
  # ------------------------------------------------------------
  if (isTRUE(return_ci)) {

    nm <- names(A_list)
    if (is.null(nm)) nm <- paste0("term", seq_len(R))

    fit_center <- vector("list", R)
    ci_pointwise <- vector("list", R)
    ci_joint     <- vector("list", R)
    df_list      <- vector("list", R)

    for (r in seq_len(R)) {
      Ar <- as.matrix(A_list[[r]])

      # center is beta_center
      fit_r <- as.numeric(Ar %*% beta_center)
      fit_center[[r]] <- fit_r

      se_r <- sehat[[r]]

      c_pt <- as.numeric(crit$pointwise_crit_adj[[r]])
      c_j  <- as.numeric(crit$joint_adj[[r]])

      lo_pt <- fit_r - c_pt * se_r
      hi_pt <- fit_r + c_pt * se_r

      lo_j  <- fit_r - c_j * se_r
      hi_j  <- fit_r + c_j * se_r

      ci_pointwise[[r]] <- cbind(lower = lo_pt, upper = hi_pt)
      ci_joint[[r]]     <- cbind(lower = lo_j,  upper = hi_j)

      df_list[[r]] <- data.frame(
        term = nm[r],
        s = seq_len(length(fit_r)),
        fit = fit_r,
        se  = se_r,
        lower_pt = lo_pt,
        upper_pt = hi_pt,
        lower_joint = lo_j,
        upper_joint = hi_j,
        edf = edf_r[r],
        df  = df_r[r],
        adj = adj_r[r],
        pointwise_crit = c_pt,
        joint_crit = c_j
      )
    }

    names(fit_center)  <- nm
    names(ci_pointwise) <- nm
    names(ci_joint)     <- nm
    names(df_list)      <- nm

    crit$ci <- list(
      fit = fit_center,         # centered at beta_center
      se  = sehat,
      ci_pointwise = ci_pointwise,
      ci_joint     = ci_joint,
      df = df_list,
      info = list(alpha = alpha, t_adjust = t_adjust,
                  centered_at = "beta_center",
                  pivot_center = "beta_hat")
    )

    # also store pivot center fits for debugging if you want
    crit$ci$fit_pivot <- fhat_pivot
  }

  # diagnostics
  if (isTRUE(return_T)) {
    attr(crit, "T_list") <- T_list
    attr(crit, "M_mat")  <- M_mat
  }
  if (isTRUE(return_bias)) {
    crit$bias_sum <- f_sum
  }

  crit
}

#' @keywords internal
#' @noRd
make_ci_newdata <- function(gamobj,
                            fn_domain = "yindex.vec",
                            grid = NULL,
                            grid.size = NULL) {

  # Get a clean model frame with predictors only
  tt <- stats::delete.response(stats::terms(gamobj))
  mf <- stats::model.frame(tt, data = gamobj$model, na.action = stats::na.pass)

  # Force plain data.frame (removes model.frame class/attrs that upset predict.gam)
  mf <- as.data.frame(mf, stringsAsFactors = TRUE)

  if(!(fn_domain %in% names(mf))) {
    stop("fn_domain '", fn_domain, "' not found among predictors: ",
         paste(names(mf), collapse = ", "))
  }

  # Choose grid
  xobs <- mf[[fn_domain]]
  if(is.null(grid)) {
    grid <- sort(unique(xobs))
    if(!is.null(grid.size) && length(grid) > grid.size) {
      grid <- grid[round(seq(1, length(grid), length.out = grid.size))]
    }
  } else {
    grid <- sort(unique(grid))
  }

  L <- length(grid)

  # Make a single typical row (fresh data.frame row)
  one <- mf[1, , drop = FALSE]
  for(nm in names(one)) {
    if(nm == fn_domain) next
    v <- mf[[nm]]

    if(is.numeric(v)) {
      one[[nm]] <- mean(v, na.rm = TRUE)
    } else if(is.factor(v)) {
      one[[nm]] <- factor(levels(v)[1], levels = levels(v))
    } else if(is.logical(v)) {
      one[[nm]] <- FALSE
    } else if(is.character(v)) {
      one[[nm]] <- v[which(nzchar(v))[1]]
      if(is.na(one[[nm]])) one[[nm]] <- v[1]
    } else {
      one[[nm]] <- v[1]
    }
  }

  # Replicate over grid
  dd <- one[rep(1, L), , drop = FALSE]
  dd[[fn_domain]] <- grid

  # Ensure plain data.frame (again)
  dd <- as.data.frame(dd, stringsAsFactors = TRUE)
  dd
}

#-------------------------------------------------------
# helper function for bootstrap-based CIs
# Build a list of linear maps A_r so that f_r(newdata) = A_r %*% beta
# Works generically for mgcv gam objects.
#
# Recommended "target": type="iterms" if you want the same objects as your current eval code.
# But we build them from lpmatrix to avoid predict() inside the bootstrap.
# Robust, CRAN-friendly: use lpmatrix attributes, not summary.gam()
# Build A matrices for (Intercept + each parametric column + each smooth term)
# using only stable gam object internals: nsdf and smooth[[k]]$first.para/last.para.
#
# Returns A_list where each A is nrow(newdata) x p (p = length(beta)).
# Then f_r(newdata) = A_list[[r]] %*% beta.
# Build A matrices that map basis-space beta -> functional coefficient curves over the
# functional domain, in a generic (CRAN-friendly) way for mgcv/pffr-style models.
#
# Returns a named list A_list of length:
#   1 + (# unique by-variables among smooths involving fn_domain)
# with:
#   A_list[["fn_intercept"]] : parametric intercept (if present) + all fn_domain smooths with no by
#   A_list[[byvar]]          : sum of all fn_domain smooths with by == byvar
#
# Each A is (nrow(newdata) x p_basis), where p_basis = length(coef(gamobj)).
#
# Notes:
# - newdata is coerced to a plain data.frame to avoid mgcv's "newdata is a model.frame" error.
# - Smooths lacking first.para/last.para are skipped (warned once).
#
#' @keywords internal
#' @noRd
build_A_functional_by <- function(gamobj,
                                  newdata,
                                  fn_domain = "yindex.vec",
                                  include_param_intercept = TRUE,
                                  strict = FALSE) {

  # plain data.frame
  newdata <- as.data.frame(newdata, stringsAsFactors = TRUE)
  class(newdata) <- "data.frame"
  attr(newdata, "terms") <- NULL

  req <- all.vars(stats::delete.response(stats::terms(gamobj)))
  miss <- setdiff(req, names(newdata))
  if(length(miss) > 0) stop("newdata missing: ", paste(miss, collapse=", "))

  X <- mgcv::predict.gam(gamobj, newdata = newdata, type = "lpmatrix")
  p_basis <- ncol(X)

  sm_list <- gamobj$smooth
  if(is.null(sm_list) || length(sm_list) == 0) stop("No smooths in model.")

  labels <- vapply(sm_list, \(s) s$label, character(1))

  # smooths that involve fn_domain (via sm$term)
  is_fn <- vapply(sm_list, function(sm) fn_domain %in% sm$term, logical(1))
  fn_idx <- which(is_fn)
  if(length(fn_idx) == 0) stop("No smooth terms include fn_domain = '", fn_domain, "'.")

  # intercept column
  int_col <- NULL
  if(include_param_intercept) {
    cn <- colnames(X)
    if(!is.null(cn) && "(Intercept)" %in% cn) int_col <- which(cn=="(Intercept)")[1]
    else if(strict) stop("No (Intercept) column found.")
    else int_col <- 1L
  }

  skipped <- character(0)
  add_smooth <- function(A, sm) {
    fp <- sm$first.para; lp <- sm$last.para
    if(is.null(fp) || is.null(lp) || length(fp)==0 || length(lp)==0) {
      skipped <<- c(skipped, sm$label)
      return(A)
    }
    cols <- fp:lp
    cols <- cols[cols>=1 & cols<=p_basis]
    if(length(cols)==0) { skipped <<- c(skipped, sm$label); return(A) }
    A[, cols] <- A[, cols] + X[, cols, drop=FALSE]
    A
  }

  # baseline fn smooths: fn_idx whose label does NOT contain ":" (heuristic)
  base_idx <- fn_idx[!grepl(":", labels[fn_idx], fixed = TRUE)]

  # by-variable names inferred from labels after ":" (e.g., "s(yindex.vec):X1")
  by_names <- unique(sub(".*:", "", labels[fn_idx][grepl(":", labels[fn_idx], fixed=TRUE)]))
  by_names <- sort(by_names)

  # Build A_list
  A_list <- list()

  # fn_intercept = intercept + baseline smooths
  A0 <- matrix(0, nrow=nrow(X), ncol=p_basis)
  if(!is.null(int_col)) A0[, int_col] <- A0[, int_col] + X[, int_col]
  if(length(base_idx)>0) for(k in base_idx) A0 <- add_smooth(A0, sm_list[[k]])
  A_list[[1]] <- A0
  names(A_list)[1] <- "fn_intercept"

  # by curves
  for(bv in by_names){
    idx <- fn_idx[grepl(paste0(":", bv), labels[fn_idx], fixed = TRUE)]
    Aj <- matrix(0, nrow=nrow(X), ncol=p_basis)
    for(k in idx) Aj <- add_smooth(Aj, sm_list[[k]])
    A_list[[length(A_list)+1]] <- Aj
    names(A_list)[length(A_list)] <- bv
  }

  if(length(skipped)>0 && strict) stop("Skipped smooths: ", paste(unique(skipped), collapse="; "))
  if(length(skipped)>0 && !strict) warning("Skipped smooths: ", paste(unique(skipped), collapse="; "))

  A_list
}


# -----------------------------------------------------------------------------
# build CI objects (same structure as wild_studentized_boot_crit$ci)
# but using:
#   - pointwise Wald criticals (z, optionally EDF/t-adjusted),
#   - joint criticals supplied as qn (from joint.qn / joint.np / joint.basis.np).
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.ci_build_wald <- function(
    glmfit,
    wi_list,
    penalty_diag,
    beta_center,
    fn_domain = "yindex.vec",
    alpha = 0.05,
    qn_joint = NULL,                      # numeric vector (term-wise) OR NULL
    qn_label = NULL,                      # for info/debug
    # match wild behavior
    t_adjust = c("both", "joint", "none", "pointwise"),
    df_floor = 1,
    se_floor_rel = 1e-6,
    A_col_tol = 0
) {

  t_adjust <- match.arg(t_adjust)

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0,1).")
  }

  if (is.null(glmfit$Vp)) stop("glmfit$Vp (variance matrix) is required for Wald CIs.")
  Sigma_hat <- as.matrix(glmfit$Vp)
  p <- ncol(Sigma_hat)
  if (!all(dim(Sigma_hat) == c(p, p))) stop("glmfit$Vp must be p x p.")

  beta_center <- as.numeric(beta_center)
  if (length(beta_center) != p) stop("beta_center must have length = ncol(glmfit$Vp).")

  if (!is.list(wi_list) || length(wi_list) < 2L) {
    stop("wi_list must be a list of length >= 2 (clusters) to compute EDF-based adjustments.")
  }

  P <- as.matrix(penalty_diag)
  if (!all(dim(P) == c(p, p))) stop("penalty_diag must be p x p.")
  P <- 0.5 * (P + t(P))

  # --- build ci_newdata + A_list exactly as in wild branch ---
  dd_ci <- make_ci_newdata(glmfit, fn_domain = fn_domain, grid.size = NULL)
  dd_ci <- droplevels(dd_ci)

  # Force numeric by-variables to 1 so A_list gives coefficient functions
  by_vars <- unique(na.omit(sapply(glmfit$smooth, function(s) s$by)))
  by_vars <- intersect(by_vars, names(dd_ci))
  for (bv in by_vars) {
    if (is.numeric(dd_ci[[bv]])) dd_ci[[bv]] <- 1
  }

  A_list <- build_A_functional_by(glmfit, newdata = dd_ci, fn_domain = fn_domain)

  nr <- vapply(A_list, nrow, integer(1))
  if (length(unique(nr)) != 1) stop("A_list matrices have inconsistent nrow().")
  m <- nr[1]
  R <- length(A_list)

  nm <- names(A_list)
  if (is.null(nm)) nm <- paste0("term", seq_len(R))

  # --- EDF per A_list element (prefer your existing helper) ---
  edf_r <- NULL
  edf_obj <- NULL

  if (exists("fgee_effective_df", mode = "function")) {
    edf_obj <- tryCatch(
      fgee_effective_df(glmfit, A_list = A_list, wi_list = wi_list, penalty_mat = P),
      error = function(e) NULL
    )
  }

  if (!is.null(edf_obj) && !is.null(edf_obj$edf_by_term)) {
    edf_r <- as.numeric(edf_obj$edf_by_term)
    if (length(edf_r) != R) edf_r <- rep(NA_real_, R)
  }

  # fallback EDF approximation if needed (keeps CRAN robustness)
  if (is.null(edf_r) || any(!is.finite(edf_r))) {
    n_clust <- length(wi_list)
    Wbar <- Reduce(`+`, lapply(wi_list, function(Wi) as.matrix(Wi))) / n_clust
    Wbar <- 0.5 * (Wbar + t(Wbar))

    # Use solve_pd if available; else base solve
    Winv <- if (exists("solve_pd", mode = "function")) solve_pd(Wbar + P) else solve(Wbar + P)
    S_mat <- Winv %*% Wbar
    dS <- diag(S_mat)

    edf_r2 <- numeric(R)
    for (r in seq_len(R)) {
      Ar <- as.matrix(A_list[[r]])
      active <- which(colSums(abs(Ar)) > A_col_tol)
      if (!length(active)) active <- seq_len(p)
      edf_r2[r] <- sum(dS[active])
    }
    edf_r <- edf_r2
  }

  n_clust <- length(wi_list)
  df_r <- pmax(df_floor, n_clust - edf_r)
  adj_r <- stats::qt(1 - alpha/2, df = df_r) / stats::qnorm(1 - alpha/2)

  # --- align qn to A_list terms ---
  .align_qn_to_A_list <- function(qn, glmfit, A_list, A_col_tol = 0) {
    R <- length(A_list)
    if (is.null(qn)) return(rep(NA_real_, R))

    qn <- as.numeric(qn)
    if (length(qn) == 1L) return(rep(qn, R))
    if (length(qn) == R) return(qn)

    # If qn is per-smooth, map conservatively by coefficient overlap
    S <- glmfit$smooth
    if (!is.null(S) && length(qn) == length(S)) {

      smooth_cols <- lapply(S, function(sm) {
        fp <- sm$first.para; lp <- sm$last.para
        if (is.null(fp) || is.null(lp) || length(fp) == 0L || length(lp) == 0L) {
          integer(0)
        } else {
          seq.int(fp, lp)
        }
      })

      out <- numeric(R)
      for (r in seq_len(R)) {
        Ar <- as.matrix(A_list[[r]])
        active <- which(colSums(abs(Ar)) > A_col_tol)
        if (!length(active)) active <- seq_len(ncol(Ar))

        hits <- which(vapply(smooth_cols, function(cols) {
          length(cols) > 0L && length(intersect(cols, active)) > 0L
        }, logical(1)))

        if (!length(hits)) hits <- 1L
        out[r] <- max(qn[hits], na.rm = TRUE)
      }

      warning(
        "qn length (", length(qn), ") != length(A_list) (", R, "). ",
        "Mapped qn to A_list by smooth-coefficient overlap (conservative)."
      )
      return(out)
    }

    stop(
      "qn length (", length(qn), ") does not match length(A_list) (", R, ") ",
      "and cannot be mapped (length(glmfit$smooth)=", length(glmfit$smooth), ")."
    )
  }

  joint_crit <- .align_qn_to_A_list(qn_joint, glmfit, A_list, A_col_tol = A_col_tol)

  # --- pointwise Wald critical (z); joint uses qn ---
  pointwise_crit <- rep(stats::qnorm(1 - alpha/2), R)

  # adjusted criticals (match wild's t_adjust semantics)
  pointwise_crit_adj <- pointwise_crit
  joint_crit_adj <- joint_crit

  if (t_adjust %in% c("pointwise", "both")) pointwise_crit_adj <- pointwise_crit_adj * adj_r
  if (t_adjust %in% c("joint", "both"))     joint_crit_adj <- joint_crit_adj * adj_r

  # --- build fit/se + CI bands ---
  fit_list <- vector("list", R)
  se_list  <- vector("list", R)
  ci_pointwise <- vector("list", R)
  ci_joint     <- vector("list", R)
  df_list      <- vector("list", R)

  for (r in seq_len(R)) {
    Ar <- as.matrix(A_list[[r]])
    if (ncol(Ar) != p) stop("A_list[[", r, "]] must have ncol = p.")

    fit_r <- as.numeric(Ar %*% beta_center)

    # diag(Ar Sigma Ar^T) without forming Ar Sigma Ar^T:
    tmp <- Ar %*% Sigma_hat                 # m x p
    v_r <- rowSums(tmp * Ar)                # length m
    se_r <- sqrt(pmax(v_r, 0))

    # floor
    se_r <- pmax(se_r, se_floor_rel * se_r + .Machine$double.eps)

    fit_list[[r]] <- fit_r
    se_list[[r]]  <- se_r

    c_pt <- as.numeric(pointwise_crit_adj[r])
    lo_pt <- fit_r - c_pt * se_r
    hi_pt <- fit_r + c_pt * se_r
    ci_pointwise[[r]] <- cbind(lower = lo_pt, upper = hi_pt)

    c_j <- as.numeric(joint_crit_adj[r])
    if (is.finite(c_j)) {
      lo_j <- fit_r - c_j * se_r
      hi_j <- fit_r + c_j * se_r
      ci_joint[[r]] <- cbind(lower = lo_j, upper = hi_j)
    } else {
      ci_joint[[r]] <- cbind(lower = rep(NA_real_, m), upper = rep(NA_real_, m))
    }

    df_list[[r]] <- data.frame(
      term = nm[r],
      s = seq_len(length(fit_r)),
      fit = fit_r,
      se  = se_r,
      lower_pt = lo_pt,
      upper_pt = hi_pt,
      lower_joint = ci_joint[[r]][, 1],
      upper_joint = ci_joint[[r]][, 2],
      edf = edf_r[r],
      df  = df_r[r],
      adj = adj_r[r],
      pointwise_crit = pointwise_crit_adj[r],
      joint_crit = joint_crit_adj[r]
    )
  }

  names(fit_list)     <- nm
  names(se_list)      <- nm
  names(ci_pointwise) <- nm
  names(ci_joint)     <- nm
  names(df_list)      <- nm

  # small infix helper
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # --- build crit object in the same “shape” as wild_studentized_boot_crit ---
  crit <- list(
    joint = as.list(as.numeric(joint_crit)),
    pt_lo = lapply(seq_len(R), function(r) rep(-pointwise_crit[r], m)),
    pt_hi = lapply(seq_len(R), function(r) rep( pointwise_crit[r], m)),
    pointwise_crit = as.list(as.numeric(pointwise_crit)),
    info = list(
      studentize = "none",
      B = NA_integer_,
      alpha = alpha,
      method = "wald",
      qn_label = qn_label %||% NA_character_,
      t_adjust = t_adjust,
      edf_r = edf_r,
      df_r  = df_r,
      adj_r = adj_r,
      beta_centered_at = beta_center
    )
  )

  crit$joint_adj <- as.list(as.numeric(joint_crit_adj))
  crit$pointwise_crit_adj <- as.list(as.numeric(pointwise_crit_adj))

  crit$ci <- list(
    fit = fit_list,
    se  = se_list,
    ci_pointwise = ci_pointwise,
    ci_joint     = ci_joint,
    df = df_list,
    info = list(
      alpha = alpha,
      t_adjust = t_adjust,
      centered_at = "beta_center",
      pivot_center = "none"
    )
  )

  # keep this slot for “same elements as wild”
  crit$ci$fit_pivot <- fit_list

  list(
    crit = crit,
    dd_ci = dd_ci,
    A_list = A_list
  )
}


