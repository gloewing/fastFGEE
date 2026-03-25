# --------------------------------------------------------------------------
# One-time fold prep for fast loss:
#  - convert grp to integer ids (required for C++)
#  - ensure each grp id appears in a single contiguous block (required for your C++)
#  - does NOT change the CV objective (order-invariant)
# --------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.prep_fold_holdout_data <- function(XX, YY, GRP) {
  K <- length(XX)
  stopifnot(length(YY) == K, length(GRP) == K)

  for (kk in seq_len(K)) {
    y <- as.numeric(YY[[kk]])
    X <- XX[[kk]]
    g <- GRP[[kk]]

    # integer group ids for C++ (safe for R too)
    g_int <- as.integer(factor(g))

    # Check contiguity: if any group id repeats in rle(values), it is NOT contiguous
    rv <- rle(g_int)$values
    not_contiguous <- any(duplicated(rv))

    if (not_contiguous) {
      o <- order(g_int)
      g_int <- g_int[o]
      y     <- y[o]
      X     <- X[o, , drop = FALSE]
    }

    GRP[[kk]] <- g_int
    YY[[kk]]  <- y
    XX[[kk]]  <- X
  }

  list(XX = XX, YY = YY, GRP = GRP)
}


# cross-validation
#' @keywords internal
#' @noRd
fun.gee1step.cv <- function(w, d, grid, data,
                            namesd,
                            cname_,
                            fid_ = NULL,          # <-- NEW: functional-observation id (often time or curve_id)
                            fit.initial,
                            cv = FALSE,
                            B = 1000,
                            K = 10,
                            folds.list = NULL,
                            sets = NULL,
                            seed = 1,
                            refine_eps = 1e-6,
                            loss = c("nll", "brier"),
                            clip_prob = 1e-6,
                            rule = c("min", "1se"),
                            se_mult = 1,
                            zero_floor = 1e-8,
                            exact = FALSE,
                            d_exact = NULL,
                            eval_prop = 1) {

  # -----------------------
  # args + basic setup
  # -----------------------
  loss <- match.arg(loss)
  rule <- select_rule <- match.arg(rule)

  if (!is.numeric(clip_prob) || length(clip_prob) != 1L) clip_prob <- 1e-6
  clip_prob <- max(min(clip_prob, 0.499999), 0)

  if (!is.numeric(se_mult) || length(se_mult) != 1L || se_mult <= 0) se_mult <- 1
  if (!is.numeric(refine_eps) || length(refine_eps) != 1L || refine_eps <= 0) refine_eps <- 1e-6

  if (!is.numeric(eval_prop) || length(eval_prop) != 1L || !is.finite(eval_prop)) eval_prop <- 1
  if (eval_prop <= 0 || eval_prop > 1) stop("eval_prop must be in (0, 1]. Got ", eval_prop, ".")

  if (eval_prop < 1 && loss != "nll") {
    warning("eval_prop < 1 is only applied when loss='nll'. Using full evaluation data for loss='", loss, "'.")
  }

  data <- data.table::copy(data.table::as.data.table(data))

  # Pull family/link
  family_raw  <- fit.initial$family$family
  family.dist <- tolower(family_raw)
  link.fn     <- tolower(fit.initial$family$link)

  # ---- evaluation subsample flag ----
  need_eval_subsample <- (loss == "nll" && eval_prop < 1)

  # Validate fid_ if needed
  fid_vec <- NULL
  if (need_eval_subsample) {
    if (is.null(fid_)) stop("eval_prop < 1 with loss='nll' requires fid_ (functional-observation id column name).")
    fid_ <- as.character(fid_)
    if (length(fid_) != 1L) stop("fid_ must be a single column name string.")
    if (!(fid_ %in% names(data))) stop("fid_ column not found in data: ", fid_)
    fid_vec <- data[[fid_]]
  }

  # solver helper
  # .solve_chol <- function(a, b) {
  #   if (requireNamespace("sanic", quietly = TRUE)) sanic::solve_chol(a = a, b = b) else solve(a, b)
  # }

  # ---- Exact Gaussian guard + choose which d to use ----
  if (exact) {
    if (!(tolower(fit.initial$family$family) == "gaussian" && link.fn == "identity")) {
      stop("estimator='exact' is only supported for Gaussian identity. Got family='",
           family.dist, "', link='", link.fn, "'.")
    } else {
      estimator <- "exact"
    }
  } else {
    estimator <- "one-step"
  }

  if (estimator == "exact") {
    if (!is.null(d_exact)) {
      d_use <- d_exact
    } else {
      warning("estimator='exact' but d_exact is NULL: using argument 'd'. ",
              "Make sure you passed EXACT Gaussian D contributions as d_exact.")
      d_use <- d
    }
  } else {
    d_use <- d
  }

  # ---- cluster vector + alignment ----
  clust.vec0 <- as.vector(data[, ..cname_])[[1]]
  clust.vec  <- as.character(clust.vec0)
  clusters   <- unique(clust.vec)
  n <- length(clusters)
  K <- min(c(n, K))
  N <- nrow(data)

  # Reorder w/d lists to match cluster order in data, if possible
  if (!is.null(names(w)) && all(clusters %in% names(w))) w <- w[clusters]
  if (!is.null(names(d_use)) && all(clusters %in% names(d_use))) d_use <- d_use[clusters]

  beta0 <- fit.initial$coefficients

  # --------------------------------------------------------------------------
  # Helper: safe grid refinement (prevents zero collapse)
  # --------------------------------------------------------------------------
  .refine_grid_safe <- function(lambda_star, base, grid_now, zero_floor = 1e-8) {
    lam <- as.numeric(lambda_star)
    dlam <- length(lam)

    for (j in seq_len(dlam)) {
      if (!is.finite(lam[j])) stop("Non-finite lambda_star encountered at dim ", j)
      if (lam[j] < zero_floor) {
        colj <- grid_now[, j]
        minpos <- suppressWarnings(min(colj[colj > 0], na.rm = TRUE))
        if (!is.finite(minpos)) minpos <- zero_floor
        lam[j] <- max(minpos, zero_floor)
      }
    }

    if (is.matrix(base)) {
      if (ncol(base) == 1L) {
        mults <- replicate(dlam, as.numeric(base[, 1]), simplify = FALSE)
      } else {
        base_df <- as.data.frame(base)
        if (ncol(base_df) != dlam) stop("Next-stage multiplier matrix must have ncol = length(lambda_star).")
        mults <- lapply(base_df, as.numeric)
      }
    } else if (is.data.frame(base)) {
      if (ncol(base) != dlam) stop("Next-stage multiplier data.frame must have ncol = length(lambda_star).")
      mults <- lapply(base, as.numeric)
    } else if (is.list(base)) {
      if (length(base) == 1L && dlam > 1L) {
        mults <- replicate(dlam, as.numeric(base[[1]]), simplify = FALSE)
      } else {
        if (length(base) != dlam) stop("Next-stage multiplier list must have length = length(lambda_star).")
        mults <- lapply(base, as.numeric)
      }
    } else {
      mults <- replicate(dlam, as.numeric(base), simplify = FALSE)
    }

    out <- expand.grid(lapply(seq_len(dlam), function(j) lam[j] * mults[[j]]))
    as.matrix(out)
  }

  # --------------------------------------------------------------------------
  # Helper: cluster-averaged loss from contribution matrix C (nobs x G)
  # --------------------------------------------------------------------------
  .cluster_average_cols <- function(Cmat, grp_int) {
    S <- rowsum(Cmat, grp_int, reorder = FALSE)
    cnt <- rowsum(rep.int(1, nrow(Cmat)), grp_int, reorder = FALSE)
    if (is.matrix(cnt)) cnt <- cnt[, 1]
    M <- sweep(S, 1, cnt, "/")
    colMeans(M)
  }

  # --------------------------------------------------------------------------
  # Link inverse for CV (matrix-safe)
  # --------------------------------------------------------------------------
  f_link_cv <- gee_family_fns(
    family    = fit.initial$family,
    link      = fit.initial$family$link,
    clamp_eps = 0
  )

  linkinv_fn <- function(eta) {
    out <- f_link_cv$linkinv(eta)
    if (!is.null(dim(eta)) && is.null(dim(out))) dim(out) <- dim(eta)
    out
  }

  # --------------------------------------------------------------------------
  # Optional C++ acceleration (Suggests-style)
  # --------------------------------------------------------------------------
  fastkfold_cpp <- get0("fastkfold_loss_all_folds_cpp", mode = "function", inherits = TRUE)
  foldloss_cpp  <- get0("fold_loss_from_eta_cluster_cpp", mode = "function", inherits = TRUE)

  has_cpp_pkgs <- requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)

  # --------------------------------------------------------------------------
  # ROBUST C++ and R parameter setup
  # --------------------------------------------------------------------------
  fam_info <- get_family_info(fit.initial)
  family_cpp <- fam_info$family_cpp
  dispersion_cpp <- fam_info$dispersion_cpp
  link_cpp <- fam_info$link_cpp

  # Check for C++ package availability and function support
  fastkfold_cpp <- get0("fastkfold_loss_all_folds_cpp", mode = "function", inherits = TRUE)
  foldloss_cpp  <- get0("fold_loss_from_eta_cluster_cpp", mode = "function", inherits = TRUE)
  has_cpp_pkgs <- requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)

  link_supported_cpp <- link_cpp %in% c("identity", "log", "logit", "probit", "cloglog", "inverse")

  loss_supported_cpp <- !is.null(family_cpp) && (
    (family_cpp %in% c("gaussian", "binomial", "quasibinomial", "poisson", "quasipoisson", "gamma")) ||
      (family_cpp %in% c("negbinomial", "beta") && loss == "nll")
  )

  use_cpp_fastkfold_loss <- has_cpp_pkgs && is.function(fastkfold_cpp) && link_supported_cpp && loss_supported_cpp
  use_cpp_foldloss       <- has_cpp_pkgs && is.function(foldloss_cpp)  && link_supported_cpp && loss_supported_cpp

  # --------------------------------------------------------------------------
  # One-time fold prep for C++:
  #  - grp must be integer and contiguous-by-grp for your run-length averaging
  # --------------------------------------------------------------------------
  .prep_folds_for_cpp_once <- function(XX, YY, GRP_int) {
    K0 <- length(XX)
    XX2 <- vector("list", K0)
    YY2 <- vector("list", K0)
    GRP2 <- vector("list", K0)

    for (kk in seq_len(K0)) {
      X <- XX[[kk]]
      y <- as.numeric(YY[[kk]])
      g_int <- GRP_int[[kk]]

      rv <- rle(g_int)$values
      need_reorder <- any(duplicated(rv))

      if (need_reorder) {
        o <- order(g_int)
        X <- X[o, , drop = FALSE]
        y <- y[o]
        g_int <- g_int[o]
      }

      XX2[[kk]] <- X
      YY2[[kk]] <- y
      GRP2[[kk]] <- g_int
    }

    list(XX = XX2, YY = YY2, GRP = GRP2)
  }

  # --------------------------------------------------------------------------
  # R loss (supports expanded families/links)
  # --------------------------------------------------------------------------
  fam_key_simple <- tolower(gsub("[^a-z0-9]", "", family.dist))
  tiny_mu <- 1e-12

  # Use dispersion_cpp for both NB theta and Beta phi
  is_nb <- (grepl("negativebinomial", fam_key_simple, fixed = TRUE) ||
              fam_key_simple %in% c("nb", "negbinomial", "negativebinomial"))

  is_beta <- fam_key_simple %in% c("beta","betar","betaregression") ||
    grepl("betareg", fam_key_simple, fixed = TRUE)

  # For NB: dispersion_cpp holds theta
  # For Beta: dispersion_cpp holds phi
  # For other families: dispersion_cpp is 1.0 (default, not used)

  .fold_loss_from_eta_cluster_R <- function(eta_mat, y, grp_int) {
    if (nrow(eta_mat) != length(y)) stop("eta_mat rows must match length(y).")
    y <- as.numeric(y)

    if (fam_key_simple %in% c("gaussian", "normal")) {
      mu <- linkinv_fn(eta_mat)
      C  <- (mu - y)^2
      return(.cluster_average_cols(C, grp_int))
    }

    if (fam_key_simple %in% c("binomial", "quasibinomial")) {
      p <- linkinv_fn(eta_mat)
      if (clip_prob > 0) p <- pmin(pmax(p, clip_prob), 1 - clip_prob)

      C <- if (loss == "brier") (p - y)^2 else -(y * log(p) + (1 - y) * log(1 - p))
      return(.cluster_average_cols(C, grp_int))
    }

    if (fam_key_simple %in% c("poisson", "quasipoisson")) {
      mu <- linkinv_fn(eta_mat)
      mu_nll <- pmax(mu, tiny_mu)

      C <- if (loss == "brier") (mu - y)^2 else (mu_nll - y * log(mu_nll))
      return(.cluster_average_cols(C, grp_int))
    }

    if (fam_key_simple %in% c("gamma", "quasigamma")) {
      mu <- linkinv_fn(eta_mat)
      mu_nll <- pmax(mu, tiny_mu)

      C <- if (loss == "brier") (mu - y)^2 else (log(mu_nll) + (y / mu_nll) - 1)
      return(.cluster_average_cols(C, grp_int))
    }

    if (fam_key_simple %in% c("inversegaussian", "inversegaussianfamily")) {
      mu <- linkinv_fn(eta_mat)
      mu <- pmax(mu, tiny_mu)
      y_pos <- pmax(y, tiny_mu)

      C <- if (loss == "brier") (mu - y)^2 else ((y - mu)^2) / (mu^2 * y_pos)
      return(.cluster_average_cols(C, grp_int))
    }

    if (is_nb) {
      mu <- linkinv_fn(eta_mat)
      mu_nll <- pmin(pmax(mu, 1e-10), 1e10) # using this to match C++ nb implementation
      #mu_nll <- pmax(mu, tiny_mu)

      if (loss == "brier") {
        C <- (mu - y)^2
        return(.cluster_average_cols(C, grp_int))
      }

      # Use dispersion_cpp which holds theta from get_family_info()
      if (!is.finite(dispersion_cpp) || dispersion_cpp <= 0) {
        stop("NB CV loss='nll' requires theta. get_family_info() returned invalid dispersion_cpp=",
             dispersion_cpp)
      }

      th <- dispersion_cpp
      C <- (th + y) * log(th + mu_nll) - y * log(mu_nll)
      return(.cluster_average_cols(C, grp_int))
    }

    if (is_beta) {
      mu <- linkinv_fn(eta_mat)
      mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
      yb <- pmin(pmax(y, 1e-6), 1 - 1e-6)

      if (loss == "brier") {
        C <- (mu - yb)^2
        return(.cluster_average_cols(C, grp_int))
      }

      # Use dispersion_cpp which holds phi from get_family_info()
      if (!is.finite(dispersion_cpp) || dispersion_cpp <= 0) {
        stop("Beta CV loss='nll' requires precision/phi. get_family_info() returned invalid dispersion_cpp=",
             dispersion_cpp)
      }

      phi <- dispersion_cpp
      a <- mu * phi
      b <- (1 - mu) * phi
      C <- lgamma(a) + lgamma(b) - lgamma(a + b) -
        (a - 1) * log(yb) - (b - 1) * log(1 - yb)

      return(.cluster_average_cols(C, grp_int))
    }

    stop("Unsupported family.dist for CV loss: ", family.dist)
  }

  # --------------------------------------------------------------------------
  # selection helper (unchanged)
  # --------------------------------------------------------------------------
  .select_lambda_from_final_stage <- function(mse_final, se_final, traceP, rule, se_mult) {
    i_min <- which.min(mse_final)
    if (rule == "min" || length(mse_final) == 1L) {
      mm <- min(mse_final)
      idx <- which(mse_final == mm)
      i_star <- if (length(idx) > 1L) idx[which.max(traceP[idx])] else i_min
      return(list(i_min = i_min, i_star = i_star))
    }

    se_min <- se_final[i_min]
    if (!is.finite(se_min)) return(list(i_min = i_min, i_star = i_min))

    thr <- mse_final[i_min] + se_mult * se_min
    ok <- which(mse_final <= thr)
    if (length(ok) == 0L) return(list(i_min = i_min, i_star = i_min))

    max_tr <- max(traceP[ok])
    ok2 <- ok[traceP[ok] == max_tr]
    i_star <- ok2[which.min(mse_final[ok2])]
    list(i_min = i_min, i_star = i_star)
  }

  # ==========================================================================
  # NOTE: other cv modes unchanged / not shown
  # ==========================================================================
  if (isTRUE(cv) || identical(cv, "fastCV")) {
    stop("Only 'kfold' and 'fastkfold' branches were updated here.")
  }

  # penalty setup once
  pen_setup <- penalty_setup(fit.initial, unpenalized = fit.initial$nsdf)

  # --------------------------------------------------------------------------
  # Precompute row indices per cluster ONCE (fast + re-usable)
  # --------------------------------------------------------------------------
  clust_rows_map <- split(seq_len(N), clust.vec)  # names sorted, but we reorder below
  clust_rows_by_id <- clust_rows_map[clusters]    # now aligned with `clusters` order
  rm(clust_rows_map)

  # --------------------------------------------------------------------------
  # Helper: get holdout rows (train = full heldout clusters)
  #         eval  = optionally subsample functional observations (cluster,fid_)
  # --------------------------------------------------------------------------
  .holdout_rows_train <- function(sets, K) {
    lapply(seq_len(K), function(kk) unlist(clust_rows_by_id[ sets[[kk]] ], use.names = FALSE))
  }

  .holdout_rows_eval <- function(sets, K) {
    if (!need_eval_subsample) return(.holdout_rows_train(sets, K))

    lapply(seq_len(K), function(kk) {
      set.seed(seed + 907L + kk)

      idxs <- sets[[kk]]
      unlist(lapply(idxs, function(j) {
        r <- clust_rows_by_id[[j]]
        f <- fid_vec[r]
        u <- unique(f)
        m <- length(u)
        keep <- max(1L, as.integer(ceiling(eval_prop * m)))
        keep_u <- sample(u, size = keep, replace = FALSE)
        r[f %in% keep_u]
      }), use.names = FALSE)
    })
  }

  # ==========================================================================
  # kfold
  # ==========================================================================
  if (identical(cv, "kfold")) {

    message("Smoothing parameter tuning: kfold")
    YY_full <- as.matrix(data[, c("Y"), with = FALSE])
    XX_full <- as.matrix(data[, ..namesd, with = FALSE])

    set.seed(seed)

    if (is.null(sets)) {
      indices <- 1:n
      sets <- split(indices[sample.int(n)], sort(rank(indices) %% K))
      beta <- replicate(K, beta0)
    } else {
      message("fold-specific beta provided")
      beta <- folds.list
    }

    beta <- as.matrix(beta)
    p_dim <- nrow(beta)

    # FULL heldout rows (for training scaling)
    holdout_rows_train <- .holdout_rows_train(sets, K)
    n_i_vec_train <- sapply(holdout_rows_train, length)

    # EVAL heldout rows (possibly subsampled by fid_ within each heldout cluster)
    holdout_rows_eval <- .holdout_rows_eval(sets, K)

    YY <- lapply(seq_len(K), function(kk) as.numeric(YY_full[ holdout_rows_eval[[kk]], 1]))
    XX <- lapply(seq_len(K), function(kk) XX_full[ holdout_rows_eval[[kk]], , drop = FALSE])

    # grp as integer per fold (for R + C++)
    GRP_int <- lapply(seq_len(K), function(kk) as.integer(factor(clust.vec[ holdout_rows_eval[[kk]] ])))

    # one-time prep for C++ if used
    if (isTRUE(use_cpp_foldloss)) {
      tmp <- .prep_folds_for_cpp_once(XX, YY, GRP_int)
      XX <- tmp$XX; YY <- tmp$YY; GRP_int <- tmp$GRP
      rm(tmp)
    }

    # fold weights (by # held-out clusters; unchanged by eval subsampling)
    ho.len <- sapply(sets, length)
    n.up   <- n - ho.len
    fold_w <- ho.len / sum(ho.len)

    # training sums via total - holdout (over clusters)
    d_total <- Reduce("+", d_use)
    w_total <- Reduce("+", w)

    d_hold <- lapply(seq_len(K), function(kk) Reduce("+", d_use[ sets[[kk]] ]))
    w_hold <- lapply(seq_len(K), function(kk) Reduce("+", w[ sets[[kk]] ]))

    d_train <- lapply(seq_len(K), function(kk) d_total - d_hold[[kk]])
    w_train <- lapply(seq_len(K), function(kk) w_total - w_hold[[kk]])

    # grid list handling
    if (is.list(grid)) {
      grid.ls <- grid
      grid.len <- length(grid.ls)
    } else {
      grid.ls <- list(grid)
      grid.len <- 1
    }

    mse.ls <- vector("list", grid.len)
    se.ls  <- vector("list", grid.len)

    grid_final <- NULL
    mse_mat_final <- NULL
    Pen_base_final <- NULL

    for (ge in seq_len(grid.len)) {

      message(paste("k-fold cv list:", ge, "of", grid.len))

      grid_now <- as.matrix(grid.ls[[ge]])
      G <- nrow(grid_now)

      Pen_base <- lapply(seq_len(G), function(g) penalty_from_setup(pen_setup, lambda = grid_now[g, ]))

      mse_mat <- matrix(NA_real_, nrow = G, ncol = K)

      for (kk in seq_len(K)) {

        wi  <- w_train[[kk]]
        di  <- d_train[[kk]]
        nup <- n.up[[kk]]

        ni_train <- n_i_vec_train[[kk]]
        scal <- (N - ni_train) / N

        Bhat <- matrix(NA_real_, nrow = p_dim, ncol = G)

        for (g in seq_len(G)) {
          Pen <- scal * Pen_base[[g]]

          if (estimator == "one-step") {
            wwi <- wi / nup + Pen
            ddi <- di - (Pen %*% beta[, kk]) * nup
            inf_fn <- solve_pd(a = wwi, b = ddi)
            Bhat[, g] <- beta[, kk] + as.numeric(inf_fn) / nup
          } else {
            bh <- solve_pd(a = wi + Pen, b = di)
            Bhat[, g] <- as.numeric(bh)
          }
        }

        eta_mat <- XX[[kk]] %*% Bhat

        if (isTRUE(use_cpp_foldloss)) {
          mse_mat[, kk] <- as.numeric(foldloss_cpp(
            eta_mat = eta_mat,
            y = YY[[kk]],
            grp = GRP_int[[kk]],
            family = family_cpp,
            link = link_cpp,
            loss = loss,
            clip_prob = clip_prob,
            dispersion = dispersion_cpp
          ))
        } else {
          mse_mat[, kk] <- .fold_loss_from_eta_cluster_R(
            eta_mat = eta_mat,
            y = YY[[kk]],
            grp_int = GRP_int[[kk]]
          )
        }
      }

      mse <- as.numeric(mse_mat %*% fold_w)
      se  <- as.numeric(apply(mse_mat, 1, stats::sd) / sqrt(K))

      mse.ls[[ge]] <- mse
      se.ls[[ge]]  <- se

      g_min <- which.min(mse)
      lambda_min_stage <- grid_now[g_min, , drop = FALSE]

      if (ge < grid.len) {
        grid.ls[[ge + 1]] <- .refine_grid_safe(lambda_min_stage, base = grid.ls[[ge + 1]],
                                               grid_now = grid_now, zero_floor = zero_floor)
      } else {
        grid_final <- grid_now
        mse_mat_final <- mse_mat
        Pen_base_final <- Pen_base
      }
    }

    mse_final <- as.numeric(mse_mat_final %*% fold_w)
    se_final  <- as.numeric(apply(mse_mat_final, 1, stats::sd) / sqrt(K))

    traceP <- vapply(Pen_base_final, function(P) sum(diag(P)), numeric(1))
    sel <- .select_lambda_from_final_stage(mse_final, se_final, traceP, select_rule, se_mult)

    g_min <- sel$i_min
    g_final <- sel$i_star

    lambda_min  <- grid_final[g_min, , drop = FALSE]
    lambda.star <- grid_final[g_final, , drop = FALSE]
    lambda_1se  <- if (select_rule == "1se") grid_final[g_final, , drop = FALSE] else NULL

    penalty_mat_full <- penalty_from_setup(pen_setup, lambda = lambda.star)

    return(list(
      mse = mse_final,
      se  = se_final,
      mse.ls = mse.ls,
      se.ls  = se.ls,
      grid = grid_final,
      estimator = estimator,
      loss = loss,
      clip_prob = clip_prob,
      select_rule = select_rule,
      se_mult = se_mult,
      lambda.star = lambda.star,
      lambda.star.min = lambda_min,
      lambda.star.1se = lambda_1se,
      argmin = g_min,
      argfinal = g_final,
      penalty_mat = penalty_mat_full,
      eval_prop = eval_prop,
      fid_ = fid_,
      used_cpp = isTRUE(use_cpp_foldloss)
    ))
  }

  # ==========================================================================
  # fastkfold
  # ==========================================================================
  else if (identical(cv, "fastkfold")) {

    message("Smoothing parameter tuning: fastkfold")

    YY_full <- as.matrix(data[, c("Y"), with = FALSE])
    XX_full <- as.matrix(data[, ..namesd, with = FALSE])

    set.seed(seed)

    if (is.null(sets)) {
      indices <- 1:n
      sets <- split(indices[sample.int(n)], sort(rank(indices) %% K))
      beta <- replicate(K, beta0)
    } else {
      message("fold-specific beta provided")
      beta <- folds.list
    }

    beta <- as.matrix(beta)
    p_dim <- nrow(beta)

    # FULL heldout rows (for training scaling)
    holdout_rows_train <- .holdout_rows_train(sets, K)
    samp_vec_train <- sapply(holdout_rows_train, length)

    # EVAL heldout rows (possibly subsampled by fid_)
    holdout_rows_eval <- .holdout_rows_eval(sets, K)

    YY <- lapply(seq_len(K), function(kk) as.numeric(YY_full[ holdout_rows_eval[[kk]], 1]))
    XX <- lapply(seq_len(K), function(kk) XX_full[ holdout_rows_eval[[kk]], , drop = FALSE])
    GRP_int <- lapply(seq_len(K), function(kk) as.integer(factor(clust.vec[ holdout_rows_eval[[kk]] ])))

    # one-time prep for C++ if used
    if (isTRUE(use_cpp_fastkfold_loss)) {
      tmp <- .prep_folds_for_cpp_once(XX, YY, GRP_int)
      XX <- tmp$XX; YY <- tmp$YY; GRP_int <- tmp$GRP
      rm(tmp)
    }

    # fold weights (by # held-out clusters; unchanged by eval subsampling)
    ho.len <- sapply(sets, length)
    fold_w <- ho.len / sum(ho.len)

    # ---- Precompute D objects for each fold via total - holdout (fast)
    d_total <- Reduce("+", d_use)
    d_hold  <- lapply(seq_len(K), function(kk) Reduce("+", d_use[ sets[[kk]] ]))

    if (estimator == "one-step") {
      D_list <- lapply(seq_len(K), function(kk) {
        di <- d_total - d_hold[[kk]]
        di * N / (N - samp_vec_train[kk])
      })
      D <- as.matrix(do.call(cbind, D_list))
      rm(D_list)
    } else {
      Dbar_list <- lapply(seq_len(K), function(kk) {
        di <- d_total - d_hold[[kk]]
        di_adj <- di * N / (N - samp_vec_train[kk])
        di_adj / n
      })
      Dbar <- as.matrix(do.call(cbind, Dbar_list))
      rm(Dbar_list)
    }

    wi <- Reduce("+", w) / n

    # grid list handling
    if (is.list(grid)) {
      grid.ls <- grid
      grid.len <- length(grid.ls)
    } else {
      grid.ls <- list(grid)
      grid.len <- 1
    }

    mse.ls <- vector("list", grid.len)
    se.ls  <- vector("list", grid.len)

    grid_final <- NULL
    mse_mat_final <- NULL
    Pen_list_final <- NULL

    for (ge in seq_len(grid.len)) {

      message(paste("fast k-fold cv list:", ge, "of", grid.len))

      grid_now <- as.matrix(grid.ls[[ge]])
      G <- nrow(grid_now)

      Pen_list <- lapply(seq_len(G), function(g) penalty_from_setup(pen_setup, lambda = grid_now[g, ]))

      delta_arr <- array(0, dim = c(p_dim, K, G))

      for (g in seq_len(G)) {
        Pen <- Pen_list[[g]]

        if (estimator == "one-step") {
          penalty_beta <- (Pen %*% beta) * n
          B_rhs <- D - penalty_beta
          inf_mat <- solve_pd(a = wi + Pen, b = B_rhs)
          delta_arr[, , g] <- inf_mat / n
        } else {
          bh_mat <- solve_pd(a = wi + (Pen / n), b = Dbar)
          delta_arr[, , g] <- bh_mat - beta
        }
      }

      # ---------- LOSS EVAL ----------
      if (isTRUE(use_cpp_fastkfold_loss)) {

        mse_mat <- as.matrix(fastkfold_cpp(
          XX_list = XX,
          YY_list = YY,
          GRP_list = GRP_int,
          beta = beta, # beta_ = beta,
          delta_arr = delta_arr,
          family = family_cpp,
          link = link_cpp,
          loss = loss,
          clip_prob = clip_prob,
          dispersion = dispersion_cpp
        ))

      } else {

        mse_mat <- matrix(NA_real_, nrow = G, ncol = K)
        for (kk in seq_len(K)) {
          Bhat_kk <- beta[, kk] + matrix(delta_arr[, kk, ], nrow = p_dim, ncol = G)
          eta_mat <- XX[[kk]] %*% Bhat_kk

          mse_mat[, kk] <- .fold_loss_from_eta_cluster_R(
            eta_mat = eta_mat,
            y = YY[[kk]],
            grp_int = GRP_int[[kk]]
          )
        }
      }

      mse <- as.numeric(mse_mat %*% fold_w)
      se  <- as.numeric(apply(mse_mat, 1, stats::sd) / sqrt(K))

      mse.ls[[ge]] <- mse
      se.ls[[ge]]  <- se

      g_min <- which.min(mse)
      lambda_min_stage <- grid_now[g_min, , drop = FALSE]

      if (ge < grid.len) {
        grid.ls[[ge + 1]] <- .refine_grid_safe(lambda_min_stage, base = grid.ls[[ge + 1]],
                                               grid_now = grid_now, zero_floor = zero_floor)
      } else {
        grid_final <- grid_now
        mse_mat_final <- mse_mat
        Pen_list_final <- Pen_list
      }
    }

    mse_final <- as.numeric(mse_mat_final %*% fold_w)
    se_final  <- as.numeric(apply(mse_mat_final, 1, stats::sd) / sqrt(K))

    traceP <- vapply(Pen_list_final, function(P) sum(diag(P)), numeric(1))
    sel <- .select_lambda_from_final_stage(mse_final, se_final, traceP, select_rule, se_mult)

    g_min <- sel$i_min
    g_final <- sel$i_star

    lambda_min  <- grid_final[g_min, , drop = FALSE]
    lambda.star <- grid_final[g_final, , drop = FALSE]
    lambda_1se  <- if (select_rule == "1se") grid_final[g_final, , drop = FALSE] else NULL

    penalty_mat_full <- penalty_from_setup(pen_setup, lambda = lambda.star)

    return(list(
      mse = mse_final,
      se  = se_final,
      mse.ls = mse.ls,
      se.ls  = se.ls,
      grid = grid_final,
      estimator = estimator,
      loss = loss,
      clip_prob = clip_prob,
      select_rule = select_rule,
      se_mult = se_mult,
      lambda.star = lambda.star,
      lambda.star.min = lambda_min,
      lambda.star.1se = lambda_1se,
      argmin = g_min,
      argfinal = g_final,
      penalty_mat = penalty_mat_full,
      eval_prop = eval_prop,
      fid_ = fid_,
      used_cpp = isTRUE(use_cpp_fastkfold_loss)
    ))
  }

  stop("Unsupported cv mode in fun.gee1step.cv(): ", cv,
       " (supported: 'kfold', 'fastkfold').")
}

#' @keywords internal
#' @noRd
fullyItr.cv <- function(w,
                        d,
                        grid,
                        data,
                        namesd,
                        X_ = namesd,
                        cname_ = "cname_",
                        fit.initial,
                        corr_fn = "independent",
                        corr_long = "independent",
                        index_fn = "yindex.vec",
                        index_long = "time",
                        rho.smooth = FALSE,
                        ar = c("mom", "yw"),
                        clamp = 0.999,
                        # CV controls
                        cv = TRUE,
                        K = 10,
                        seed = 1,
                        sets = NULL,
                        folds.list = NULL,
                        # Iterated GEE controls (training fits inside CV)
                        beta.tol = 1e-3,
                        gee.maxiter = 5,
                        linpred_method = c("accumulate", "matrix"),
                        clip_mu = 0,
                        # penalty scaling (kept from your old code)
                        scale_lambda_by_rows = TRUE,
                        # loss controls (NEW but default matches your one-step CV)
                        loss = c("nll", "brier"),
                        clip_prob = 1e-6,
                        # misc
                        verbose = FALSE) {

  ar <- match.arg(ar)
  linpred_method <- match.arg(linpred_method)
  loss <- match.arg(loss)

  # ---- data copy (local only) ----
  dt0 <- data.table::copy(data.table::as.data.table(data))

  # canonical ordering once (safe; subset preserves order)
  if (all(c(cname_, index_long, index_fn) %in% names(dt0))) {
    data.table::setorderv(dt0, c(cname_, index_long, index_fn))
  } else {
    stop("fullyItr.cv: data must contain cname_/index_long/index_fn columns: ",
         paste(cname_, index_long, index_fn, collapse = ", "))
  }

  # ---- basic checks ----
  if (!("Y" %in% names(dt0))) stop("fullyItr.cv expects outcome column named 'Y' in long data.")
  if (!all(X_ %in% names(dt0))) stop("Some X_ columns not found in data.")
  if (length(fit.initial$coefficients) != length(X_)) {
    stop("length(fit.initial$coefficients) != length(X_).")
  }

  # align namesd vs X_ (legacy-safe)
  if (!identical(as.character(namesd), as.character(X_))) {
    if (verbose) message("fullyItr.cv: namesd != X_. Using X_ for fitting and prediction.")
    namesd <- X_
  }

  corr_fn   <- tolower(corr_fn)
  corr_long <- tolower(corr_long)

  # ---- link + family ----
  link_lc <- tolower(fit.initial$family$link)
  fam_raw <- fit.initial$family$family
  fam_key_simple <- tolower(gsub("[^a-z0-9]", "", fam_raw))

  # ------------------------------------------------------------------
  # FIX #1: define link.fn (the error you saw)
  # ------------------------------------------------------------------
  link.fn <- link_lc

  # ------------------------------------------------------------------
  # Family/nuisance extraction (same helper used by fun.gee1step.cv)
  # ------------------------------------------------------------------
  fi <- get_family_info(fit.initial)
  family_cpp     <- fi$family_cpp
  dispersion_cpp <- fi$dispersion_cpp
  link_cpp       <- fi$link_cpp

  # ------------------------------------------------------------------
  # Optional C++ loss (same style gatekeeping as fun.gee1step.cv)
  # ------------------------------------------------------------------
  foldloss_cpp <- get0("fold_loss_from_eta_cluster_cpp", mode = "function", inherits = TRUE)
  has_cpp_pkgs <- requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)

  link_supported_cpp <- link_cpp %in% c("identity", "log", "logit", "probit", "cloglog", "inverse")

  loss_supported_cpp <- !is.null(family_cpp) && (
    (family_cpp %in% c("gaussian", "binomial", "quasibinomial", "poisson", "quasipoisson", "gamma")) ||
      (family_cpp %in% c("negbinomial", "beta") && loss == "nll")
  )

  use_cpp_foldloss <- has_cpp_pkgs && is.function(foldloss_cpp) && link_supported_cpp && loss_supported_cpp

  # ------------------------------------------------------------------
  # Link inverse for R fallback loss (matrix-safe)
  # ------------------------------------------------------------------
  f_link <- gee_family_fns(
    family = fit.initial$family,
    link   = link.fn,
    clamp_eps = 0
  )

  linkinv_fn <- function(eta) {
    out <- f_link$linkinv(eta)
    if (!is.null(dim(eta)) && is.null(dim(out))) dim(out) <- dim(eta)
    out
  }

  # ------------------------------------------------------------------
  # Cluster-average helper (scalar loss for 1 eta column)
  # ------------------------------------------------------------------
  .cluster_avg_scalar <- function(Cvec, grp_int) {
    S <- rowsum(Cvec, grp_int, reorder = FALSE)
    cnt <- rowsum(rep.int(1, length(Cvec)), grp_int, reorder = FALSE)
    if (is.matrix(cnt)) cnt <- cnt[, 1]
    m <- S[, 1] / cnt
    mean(m)
  }

  # ------------------------------------------------------------------
  # R fallback loss (mirrors fun.gee1step.cv logic, incl NB/Beta)
  # ------------------------------------------------------------------
  tiny_mu <- 1e-12
  is_nb <- (grepl("negativebinomial", fam_key_simple, fixed = TRUE) ||
              fam_key_simple %in% c("nb", "negbinomial", "negativebinomial"))
  is_beta <- fam_key_simple %in% c("beta","betar","betaregression") ||
    grepl("betareg", fam_key_simple, fixed = TRUE)

  .fold_loss_scalar_R <- function(eta_vec, y, grp_int) {
    # eta_vec is numeric vector length nobs
    eta_mat <- matrix(eta_vec, ncol = 1)

    if (fam_key_simple %in% c("gaussian", "normal")) {
      mu <- as.numeric(linkinv_fn(eta_mat))
      C  <- (mu - y)^2
      return(.cluster_avg_scalar(C, grp_int))
    }

    if (fam_key_simple %in% c("binomial", "quasibinomial")) {
      p <- as.numeric(linkinv_fn(eta_mat))
      if (clip_prob > 0) p <- pmin(pmax(p, clip_prob), 1 - clip_prob)

      C <- if (loss == "brier") (p - y)^2 else -(y * log(p) + (1 - y) * log(1 - p))
      return(.cluster_avg_scalar(C, grp_int))
    }

    if (fam_key_simple %in% c("poisson", "quasipoisson")) {
      mu <- as.numeric(linkinv_fn(eta_mat))
      mu_nll <- pmax(mu, tiny_mu)
      C <- if (loss == "brier") (mu - y)^2 else (mu_nll - y * log(mu_nll))
      return(.cluster_avg_scalar(C, grp_int))
    }

    if (fam_key_simple %in% c("gamma", "quasigamma")) {
      mu <- as.numeric(linkinv_fn(eta_mat))
      mu_nll <- pmax(mu, tiny_mu)
      C <- if (loss == "brier") (mu - y)^2 else (log(mu_nll) + (y / mu_nll) - 1)
      return(.cluster_avg_scalar(C, grp_int))
    }

    if (fam_key_simple %in% c("inversegaussian", "inversegaussianfamily")) {
      mu <- as.numeric(linkinv_fn(eta_mat))
      mu <- pmax(mu, tiny_mu)
      y_pos <- pmax(y, tiny_mu)
      C <- if (loss == "brier") (mu - y)^2 else ((y - mu)^2) / (mu^2 * y_pos)
      return(.cluster_avg_scalar(C, grp_int))
    }

    if (is_nb) {
      mu <- as.numeric(linkinv_fn(eta_mat))
      mu_nll <- pmin(pmax(mu, 1e-10), 1e10)  # match your C++ NB clamp

      if (loss == "brier") {
        C <- (mu - y)^2
        return(.cluster_avg_scalar(C, grp_int))
      }

      if (!is.finite(dispersion_cpp) || dispersion_cpp <= 0) {
        stop("NB loss requires theta; get_family_info() returned invalid dispersion_cpp=", dispersion_cpp)
      }
      th <- dispersion_cpp
      C <- (th + y) * log(th + mu_nll) - y * log(mu_nll)
      return(.cluster_avg_scalar(C, grp_int))
    }

    if (is_beta) {
      mu <- as.numeric(linkinv_fn(eta_mat))
      mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
      yb <- pmin(pmax(y,  1e-6), 1 - 1e-6)

      if (loss == "brier") {
        C <- (mu - yb)^2
        return(.cluster_avg_scalar(C, grp_int))
      }

      if (!is.finite(dispersion_cpp) || dispersion_cpp <= 0) {
        stop("Beta loss requires phi; get_family_info() returned invalid dispersion_cpp=", dispersion_cpp)
      }
      phi <- dispersion_cpp
      a <- mu * phi
      b <- (1 - mu) * phi
      C <- lgamma(a) + lgamma(b) - lgamma(a + b) -
        (a - 1) * log(yb) - (b - 1) * log(1 - yb)

      return(.cluster_avg_scalar(C, grp_int))
    }

    stop("fullyItr.cv: unsupported family for loss: ", fam_raw)
  }

  # ------------------------------------------------------------------
  # Folds: match fun.gee1step.cv behavior
  #   - If sets is provided, use it
  #   - Otherwise generate the same way
  # ------------------------------------------------------------------
  clust.vec <- dt0[[cname_]]
  clusters <- unique(as.character(clust.vec))
  n_clust <- length(clusters)
  if (n_clust < 2) stop("Need at least 2 clusters for CV.")

  # treat cv as "truthy" unless explicitly FALSE
  cv_flag <- !isFALSE(cv)

  if (cv_flag) {

    if (is.null(sets)) {
      K <- min(K, n_clust)
      if (K < 2) stop("K must be >= 2 (or you need more clusters).")

      set.seed(seed)
      indices <- seq_len(n_clust)
      sets <- split(indices[sample.int(n_clust)], sort(rank(indices) %% K))
    } else {
      # use provided sets (to match one-step CV folds)
      if (!is.list(sets) || length(sets) < 2) stop("fullyItr.cv: provided 'sets' must be a list of folds (length>=2).")
      K <- length(sets)

      # allow sets to be either numeric indices or cluster IDs
      sets <- lapply(sets, function(s) {
        if (is.character(s)) {
          idx <- match(as.character(s), clusters)
          if (anyNA(idx)) stop("fullyItr.cv: some cluster IDs in sets not found in data clusters.")
          idx
        } else {
          as.integer(s)
        }
      })
    }

  } else {
    sets <- list(integer(0))
    K <- 1
  }

  # fold weights (match fun.gee1step.cv: by held-out cluster count)
  ho.len <- sapply(sets, length)
  if (sum(ho.len) > 0) fold_w <- ho.len / sum(ho.len) else fold_w <- rep(1 / K, K)

  # row indices by cluster id
  row_by_cluster <- split(seq_len(nrow(dt0)), as.character(clust.vec))

  holdout_rows <- lapply(seq_len(K), function(kk) {
    hold_ids <- as.character(clusters[sets[[kk]]])
    if (length(hold_ids) == 0) return(integer(0))
    unlist(row_by_cluster[hold_ids], use.names = FALSE)
  })

  # precompute holdout X/Y/grp once
  Y_hold <- lapply(seq_len(K), function(kk) {
    rr <- holdout_rows[[kk]]
    if (!length(rr)) return(numeric(0))
    as.numeric(dt0[rr, "Y", with = FALSE][[1]])
  })

  X_hold <- lapply(seq_len(K), function(kk) {
    rr <- holdout_rows[[kk]]
    if (!length(rr)) return(matrix(numeric(0), nrow = 0, ncol = length(namesd)))
    as.matrix(dt0[rr, ..namesd])
  })

  GRP_hold <- lapply(seq_len(K), function(kk) {
    rr <- holdout_rows[[kk]]
    if (!length(rr)) return(integer(0))
    as.character(dt0[[cname_]][rr])
  })

  # ensure integer grp + contiguity for C++ (safe even if we use R)
  tmp <- .prep_fold_holdout_data(XX = X_hold, YY = Y_hold, GRP = GRP_hold)
  # X_hold <- tmp$XX
  Y_hold <- tmp$YY
  GRP_int <- tmp$GRP
  rm(tmp)

  # ---- penalty setup ----
  pen_setup <- penalty_setup(fit.initial, unpenalized = fit.initial$nsdf)

  # ---- normalize grid into list-of-stages ----
  normalize_grid <- function(g) {
    if (is.data.frame(g)) g <- as.matrix(g)
    if (is.vector(g) && !is.list(g)) g <- matrix(g, ncol = 1)
    if (!is.matrix(g)) stop("grid stage must be a matrix/data.frame/numeric vector.")
    storage.mode(g) <- "numeric"
    if (nrow(g) < 1) stop("grid stage has zero rows.")
    g
  }

  if (is.list(grid)) {
    grid.ls <- grid
  } else {
    grid.ls <- list(grid)
  }
  grid.ls <- lapply(grid.ls, normalize_grid)
  n_stage <- length(grid.ls)

  # ---- align w/d lists with clusters ----
  if (!is.null(names(w)) && all(as.character(clusters) %in% names(w))) {
    w <- w[as.character(clusters)]
  }
  if (!is.null(names(d)) && all(as.character(clusters) %in% names(d))) {
    d <- d[as.character(clusters)]
  }
  if (length(w) != n_clust || length(d) != n_clust) {
    stop("fullyItr.cv: w/d must be lists with one element per cluster and aligned to 'clusters'.")
  }

  .eta_from_rows <- function(dt, rr, xcols, beta) {
    eta <- numeric(length(rr))
    for (j in seq_along(xcols)) {
      bj <- beta[j]
      if (!is.finite(bj) || bj == 0) next
      eta <- eta + dt[[xcols[j]]][rr] * bj
    }
    eta
  }

  # ---- one training fit for one lambda (fully iterated) ----
  fit_train_iterated <- function(dt_train, beta_start, penalty_mat, n_train_clusters) {

    beta <- as.numeric(beta_start)

    # initialize nuisance once per lambda fit
    nu <- list(
      dispersion = NULL,
      theta      = if (fi$family_cpp == "negbinomial") fi$dispersion_cpp else NULL,
      precision  = if (fi$family_cpp == "beta")       fi$dispersion_cpp else NULL,
      zi_prob    = NULL
    )

    for (it in seq_len(gee.maxiter)) {
      # message(paste("Tune itr:", it))
      dt_train <- fgee_update_working_cols_dt(
        dx = dt_train,
        namesd = X_,
        beta = beta,
        family = fit.initial$family,         # use  *real* family object
        link   = fit.initial$family$link,
        exact = FALSE,
        gaussian_v_by = index_fn,
        update_nuisance = if (it == 1L) "fixed" else "moment",
        dispersion = nu$dispersion,
        theta      = nu$theta,
        precision  = nu$precision,
        zi_prob    = nu$zi_prob,
        clamp_eps = max(1e-8, min(clip_mu, 0.499999)),
        linpred_method = linpred_method,
        eta_by = NULL,
        copy = FALSE
      )

      # persist nuisance
      nu <- attr(dt_train, "nuisance")

      cor_it <- corr.est(
        dx = dt_train,
        cname_ = cname_,
        index_fn = index_fn,
        index_long = index_long,
        corr_fn = corr_fn,
        corr_long = corr_long,
        resid_col = "resid",
        rho.smooth = rho.smooth,
        ar = ar,
        glmfit = if (isTRUE(rho.smooth)) fit.initial else NULL,
        fpca_fn = NULL,
        clamp = clamp,
        copy_dt = FALSE,
        verbose = FALSE
      )
      dt_train <- cor_it$dx
      rm(cor_it)

      wi <- W.estimate(
        dt_train,
        namesd = X_,
        cname_ = cname_,
        corr_fn = corr_fn,
        corr_long = corr_long,
        index_fn = index_fn,
        index_long = index_long,
        algo = "gschur",
        fpca_fn = cor_it$fpca,
        copy_dt = FALSE
      )

      di <- D.estimate(
        dt_train,
        namesd = X_,
        cname_ = cname_,
        corr_fn = corr_fn,
        corr_long = corr_long,
        index_fn = index_fn,
        index_long = index_long,
        algo = "gschur",
        fpca_fn = cor_it$fpca,
        copy_dt = FALSE,
        resid_col = "resid"
      )

      wi_bar <- Reduce("+", wi) / n_train_clusters
      rm(wi)
      di_sum <- Reduce("+", di)
      rm(di)

      if (it %% 2L == 0L) gc(FALSE)

      pen_vec <- as.numeric(penalty_mat %*% beta)

      step <- solve_pd(
        a = wi_bar + penalty_mat,
        b = di_sum - pen_vec * n_train_clusters
      )

      beta_new <- beta + as.numeric(step) / n_train_clusters
      if (max(abs(beta_new - beta)) <= beta.tol) {
        beta <- beta_new
        break
      }
      beta <- beta_new
    }

    beta
  }

  # ---- results per stage ----
  mse.ls <- vector("list", n_stage)
  best_lambda <- NULL

  N_full <- nrow(dt0)

  for (ge in seq_len(n_stage)) {

    grid_stage <- grid.ls[[ge]]
    grid_size <- nrow(grid_stage)

    if (verbose) message("fullyItr.cv: stage ", ge, "/", n_stage,
                         " with ", grid_size, " lambda candidates.")

    mse_mat <- matrix(NA_real_, nrow = grid_size, ncol = K)

    for (kk in seq_len(K)) {

      hold_idx <- sets[[kk]]
      train_idx <- setdiff(seq_len(n_clust), hold_idx)

      n_train_clusters <- length(train_idx)
      if (n_train_clusters < 1) stop("Fold has no training clusters (unexpected).")

      # training data view
      dt_train <- if (length(holdout_rows[[kk]]) > 0) dt0[-holdout_rows[[kk]]] else dt0

      # legacy: scale lambda by fraction of rows retained
      row_scale <- 1
      if (isTRUE(scale_lambda_by_rows)) {
        row_scale <- (N_full - length(holdout_rows[[kk]])) / N_full
      }

      # one-step warm start for FIRST lambda
      Wbar0_train <- Reduce("+", w[train_idx]) / n_train_clusters
      Dsum0_train <- Reduce("+", d[train_idx])

      beta0 <- as.numeric(fit.initial$coefficients)

      lambda1 <- as.numeric(grid_stage[1, ]) * row_scale
      P1 <- penalty_from_setup(pen_setup, lambda = lambda1)
      if (all(lambda1 == 0)) P1[,] <- 0

      pen_vec0 <- as.numeric(P1 %*% beta0)

      step1 <- solve_pd(
        a = Wbar0_train + P1,
        b = Dsum0_train - pen_vec0 * n_train_clusters
      )

      beta_start <- beta0 + as.numeric(step1) / n_train_clusters

      for (g in seq_len(grid_size)) {

        lambda_g <- as.numeric(grid_stage[g, ]) * row_scale
        P_g <- penalty_from_setup(pen_setup, lambda = lambda_g)
        if (all(lambda_g == 0)) P_g[,] <- 0

        beta_hat_g <- tryCatch(
          fit_train_iterated(dt_train = dt_train,
                             beta_start = beta_start,
                             penalty_mat = P_g,
                             n_train_clusters = n_train_clusters ),
          error = function(e) NA_real_
        )

        if (anyNA(beta_hat_g)) {
          mse_mat[g, kk] <- Inf
          next
        }

        beta_start <- beta_hat_g

        # evaluate on holdout (cluster-averaged loss, same style as fun.gee1step.cv)
        if (length(holdout_rows[[kk]]) > 0) {

          #eta <- as.numeric(X_hold[[kk]] %*% beta_hat_g)
          rr <- holdout_rows[[kk]]
          eta <- .eta_from_rows(dt0, rr, namesd, beta_hat_g)

          if (isTRUE(use_cpp_foldloss)) {

            mse_mat[g, kk] <- as.numeric(foldloss_cpp(
              eta_mat    = matrix(eta, ncol = 1),
              y          = Y_hold[[kk]],
              grp        = GRP_int[[kk]],
              family     = family_cpp,
              link       = link_cpp,
              loss       = loss,
              clip_prob  = clip_prob,
              dispersion = dispersion_cpp
            ))[1]

          } else {

            mse_mat[g, kk] <- .fold_loss_scalar_R(
              eta_vec  = eta,
              y        = Y_hold[[kk]],
              grp_int  = GRP_int[[kk]]
            )
          }

        } else {
          mse_mat[g, kk] <- NA_real_
        }

        rm(P_g, beta_hat_g, eta)
        if (g %% 5L == 0L) gc(FALSE)

      }

      rm(dt_train, Wbar0_train, Dsum0_train)
      gc(FALSE)

    }

    # fold-weighted average (match fun.gee1step.cv)
    mse <- as.numeric(mse_mat %*% fold_w)
    mse.ls[[ge]] <- list(grid = grid_stage, mse = mse, mse_mat = mse_mat)

    best_row <- which.min(mse)
    best_lambda <- as.numeric(grid_stage[best_row, , drop = TRUE])

    if (verbose) {
      message("  stage ", ge, " best lambda: ", paste(signif(best_lambda, 6), collapse = ", "),
              "  loss=", signif(mse[best_row], 6))
    }

    # refine next stage grid if there is one
    if (ge < n_stage) {

      next_base <- grid.ls[[ge + 1]]

      if (is.matrix(next_base) && ncol(next_base) == 1L) {
        mult <- as.numeric(next_base[, 1])
        refined <- expand.grid(lapply(best_lambda, function(x) x * mult))
        grid.ls[[ge + 1]] <- as.matrix(refined)
      } else if (is.matrix(next_base) && ncol(next_base) == length(best_lambda)) {
        refined <- sweep(next_base, 2, best_lambda, `*`)
        grid.ls[[ge + 1]] <- refined
      } else {
        stop("fullyItr.cv: cannot refine next grid stage. ",
             "Expected next stage to have ncol=1 or ncol=length(lambda).")
      }
    }
  }

  # final best from last stage
  final_stage <- mse.ls[[n_stage]]
  grid_final <- final_stage$grid
  mse_final  <- final_stage$mse
  best_row   <- which.min(mse_final)
  lambda.star <- as.numeric(grid_final[best_row, , drop = TRUE])

  penalty_star <- penalty_from_setup(pen_setup, lambda = lambda.star)
  if (all(lambda.star == 0)) penalty_star[,] <- 0

  list(
    mse = mse_final,
    mse.ls = lapply(mse.ls, function(z) z$mse),
    grid = lapply(mse.ls, function(z) z$grid),
    mse_mat = lapply(mse.ls, function(z) z$mse_mat),
    lambda.star = lambda.star,
    penalty_mat = penalty_star,
    bias = rep(NA_real_, length(mse_final)),
    asy.var = rep(NA_real_, length(mse_final)),
    used_cpp = isTRUE(use_cpp_foldloss)
  )
}

#' @keywords internal
#' @noRd
#-----------------------------------------------------------------------------
# New helper function for smoothing parameter tuning (with penalty setup included)
tune_smoothing_parameters <- function(
    glmfit,
    cv.grid,
    wi0,
    di0,
    di0_exact = NULL,
    use_exact_gaussian = FALSE,
    namesd,
    dx,
    tune.method = c("one-step", "fully-iterated"),
    tune.full_fn = NULL,
    folds.list = NULL,
    sets = NULL,
    preTn = TRUE,
    corr_fn = "independent",
    corr_long = "independent",
    index_fn = "yindex.vec",
    index_long = "time",
    rho.smooth = TRUE,
    X_ = NULL,
    cv = TRUE,
    eval_prop = 1,
    fid_ = "time"
) {
  tune.method <- match.arg(tune.method)

  # ------------------------------------------------------------
  # Penalty setup once per fit
  # ------------------------------------------------------------
  pen_setup <- penalty_setup(glmfit, unpenalized = glmfit$nsdf)
  sp_vec <- glmfit[["sp"]]
  if (is.list(sp_vec)) sp_vec <- sp_vec[[1]]
  sp_vec <- as.numeric(sp_vec)
  if (length(sp_vec) == 0L) sp_vec <- 0

  penalty_diag <- penalty_from_setup(pen_setup, lambda = sp_vec)
  if (all(sp_vec == 0)) {
    penalty_diag[,] <- 0
  }

  # ------------------------------------------------------------
  # Tuning (default: one-step tuning using wi0/di0)
  #   NOTE: For exact Gaussian mode, this still uses SCORE di0.
  #         You said you will update CV for exact later; this avoids
  #         changing CV behavior now.
  # ------------------------------------------------------------
  message("Smoothing Parameter Tuning")

  if (!is.null(cv.grid)) {
    # normalize cv.grid to list-of-stages format
    if (is.list(cv.grid)) {
      grid <- cv.grid
      if (length(grid) > 1L) {
        tune.ind <- TRUE
      } else {
        tune.ind <- (length(grid[[1]]) > 1L)
        if (tune.ind) grid <- list(grid, grid)
      }
    } else {
      grid <- list(cv.grid)
      tune.ind <- (length(grid[[1]]) > 1L)
    }

    if (tune.ind) {
      # baseline scaling for stage 1
      base_lambda <- sp_vec
      ok_len <- c(1L, pen_setup$n_sm, pen_setup$n_pen)

      if (!(length(base_lambda) %in% ok_len)) {
        warning("glmfit$sp has length ", length(base_lambda),
                " but penalty_setup expects length 1, #smooths=", pen_setup$n_sm,
                ", or #penalties=", pen_setup$n_pen,
                ". Using scalar baseline 1 for pre-scaling.")
        base_lambda <- 1
      }

      if (!preTn) {
        Lambda_init <- rep(1, length(base_lambda))
        grid[[1]] <- as.matrix(expand.grid(lapply(Lambda_init, function(x) x * grid[[1]])))
      } else {
        grid[[1]] <- do.call(cbind, lapply(base_lambda, function(x) x * grid[[1]]))
      }

      if (tune.method == "one-step") {
        cv.lambda <- fun.gee1step.cv(
          w = wi0,
          d = di0,
          d_exact = di0_exact,
          exact = use_exact_gaussian,
          namesd = namesd,
          grid = grid,
          data = dx,
          cname_ = "cname_",
          fit.initial = glmfit,
          folds.list = folds.list,
          loss = "nll",
          sets = sets,
          cv = "fastkfold",
          rule = "min",
          se_mult = 1,
          B = 500,
          seed = 1,
          eval_prop = eval_prop,
          fid_ = fid_
        )
      } else {
        # fully-iterated tuning hook
        if (is.null(tune.full_fn)) {
          if (exists("fullyItr.cv", mode = "function")) {
            tune.full_fn <- get("fullyItr.cv", mode = "function")
          } else {
            stop("tune.method='fully-iterated' requires tune.full_fn or a global fullyItr.cv().")
          }
        }
        cv.lambda <- tune.full_fn(
          w = wi0,
          d = di0,
          namesd = namesd,
          data = dx,
          cname_ = "cname_",
          fit.initial = glmfit,
          cv = cv,
          grid = grid,
          X_ = X_,
          corr_fn = corr_fn,
          corr_long = corr_long,
          index_fn = index_fn,
          index_long = index_long,
          rho.smooth = rho.smooth,
          sets = sets,
          folds.list = folds.list
        )
      }

      # obtain penalty matrix
      if (!is.null(cv.lambda$penalty_mat)) {
        penalty_diag <- cv.lambda$penalty_mat
      } else if (!is.null(cv.lambda$lambda.star)) {
        penalty_diag <- penalty_from_setup(pen_setup, lambda = cv.lambda$lambda.star)
      } else {
        stop("CV object must provide either $penalty_mat or $lambda.star.")
      }

      message("CV Complete")

    } else {
      # no tuning: single lambda specification
      lambda_star <- grid[[1]]
      if (is.data.frame(lambda_star)) lambda_star <- as.matrix(lambda_star)
      if (is.matrix(lambda_star)) {
        if (nrow(lambda_star) != 1L) {
          stop("Non-tuning cv.grid must represent a single lambda (vector or 1-row matrix).")
        }
        lambda_star <- as.numeric(lambda_star[1, ])
      } else {
        lambda_star <- as.numeric(lambda_star)
      }

      penalty_diag <- penalty_from_setup(pen_setup, lambda = lambda_star)
      cv.lambda <- list(
        mse = NA_real_,
        mse.ls = NULL,
        lambda.star = lambda_star,
        grid = grid,
        bias = NA_real_,
        asy.var = NA_real_,
        penalty_mat = penalty_diag
      )
    }
  } else {
    message("No penalty")
    penalty_diag <- penalty_from_setup(pen_setup, lambda = 0)
    cv.lambda <- NULL
  }

  return(list(
    penalty_diag = penalty_diag,
    cv.lambda = cv.lambda
  ))
}
