# Optimized one-step fGEE engine ----------------------------------------------
#
# The optimized engine shares one compact W/D construction pass, uses
# coefficient-space Gaussian CV sufficient statistics, and supports compact
# sandwich and wild-bootstrap inference.  The legacy engine remains available
# for comparison and unsupported workflows.

.fgee_update_mean_and_corr <- function(dx,
                                       beta,
                                       namesd,
                                       glmfit,
                                       corr_fn,
                                       corr_long,
                                       index_fn,
                                       index_long,
                                       rho.smooth,
                                       rho.pool = c("fn", "none", "long", "both"),
                                       linpred_method,
                                       clip_mu,
                                       update_nuisance = c("fixed", "moment"),
                                       nuisance = NULL) {
  update_nuisance <- match.arg(update_nuisance)
  fi <- get_family_info(glmfit)
  nuisance <- .fgee_or(nuisance, .fgee_nuisance_legacy_values(fi$nuisance))

  dx <- fgee_update_working_cols_dt(
    dx = dx,
    namesd = namesd,
    beta = beta,
    family = glmfit$family,
    link = glmfit$family$link,
    exact = FALSE,
    gaussian_v_by = index_fn,
    update_nuisance = update_nuisance,
    theta = .fgee_or(nuisance$theta,
                     if (fi$family_cpp == "negbinomial") fi$dispersion_cpp else NULL),
    precision = .fgee_or(nuisance$precision,
                         if (fi$family_cpp == "beta") fi$dispersion_cpp else NULL),
    dispersion = nuisance$dispersion,
    zi_prob = nuisance$zi_prob,
    clamp_eps = max(1e-8, min(clip_mu, 0.499999)),
    linpred_method = linpred_method,
    eta_by = NULL,
    copy = FALSE
  )

  nuisance_new <- attr(dx, "nuisance")
  cor <- corr.est(
    dx = dx,
    cname_ = "cname_",
    index_fn = index_fn,
    index_long = index_long,
    corr_fn = corr_fn,
    corr_long = corr_long,
    resid_col = "resid",
    rho.smooth = rho.smooth,
    rho.pool = rho.pool,
    ar = "mom",
    glmfit = if (isTRUE(rho.smooth)) glmfit else NULL,
    fpca_fn = NULL,
    copy_dt = FALSE
  )

  list(dx = cor$dx, rho = cor$rho, fpca = cor$fpca, nuisance = nuisance_new)
}

.fgee_one_step_update <- function(working, beta, penalty, exact = FALSE) {
  beta <- as.numeric(beta)
  P <- as.matrix(penalty)
  if (isTRUE(exact)) {
    if (is.null(working$d_exact_sum)) stop("Exact Gaussian statistics are unavailable.")
    return(as.numeric(solve_pd(working$W_sum + P, working$d_exact_sum)))
  }
  beta + as.numeric(solve_pd(
    working$W_bar + P,
    working$d_bar - as.numeric(P %*% beta)
  ))
}

.fgee_compat_scores <- function(working, exact = FALSE) {
  # Avoid duplicating the compact p x N score matrix as N separate vectors in
  # the default scores-retention path.  Full retention is the explicit
  # backwards-compatible mode for users who need the legacy list fields.
  if (!identical(working$retain, "full")) return(NULL)
  .fgee_working_dlist(working, exact = exact, required = FALSE)
}

.fgee_compat_bread <- function(working) {
  if (!identical(working$retain, "full")) return(NULL)
  .fgee_working_wlist(working, required = FALSE)
}

#' Optimized internal GEE engine
#' @keywords internal
#' @noRd
fun.gee1step.dist_itr.optimized <- function(
    orig.data,
    dx,
    formula,
    X_,
    Y_,
    namesd,
    N_clusters,
    clusters,
    glmfit,
    yindex.vec,
    bs,
    corr_fn = "independent",
    corr_long = "independent",
    index_fn = "yindex.vec",
    index_long = "time",
    var.type = "sandwich",
    cv.grid,
    boot.samps = 3000,
    rho.smooth = TRUE,
    rho.pool = c("fn", "none", "long", "both"),
    joint.CI = "wild",
    index = "yindex.vec",
    linpred_method = c("accumulate", "matrix"),
    clip_mu = 0,
    beta.tol = 1e-6,
    max.iter = 1,
    tune.method = c("one-step", "fully-iterated"),
    exact = FALSE,
    gee.fit = TRUE,
    sp.method = c(
      "auto", "fastk_staged", "fastk_grad", "sandwich_qreml",
      "qreml_fastk", "legacy"
    ),
    working.retain = c("auto", "scores", "aggregate", "full"),
    corr.solver = c("auto", "exact", "supergauss"),
    fastk.K = 10L,
    fastk.seed = 1L,
    fastk.memory = c("balanced", "speed", "lowmem"),
    fastk.start = c("qreml", "robust", "compact"),
    qreml.phi.method = c("penalized", "all", "fixed"),
    qreml.phi.fixed = NULL,
    qreml.phi.weight = 1,
    qreml.phi.clip = c(1, 8),
    keep.tuning.workspace = FALSE,
    verbose.tuning = TRUE) {

  linpred_method <- match.arg(linpred_method)
  tune.method <- match.arg(tune.method)
  working.retain <- match.arg(working.retain)
  corr.solver <- match.arg(corr.solver)
  fastk.memory <- match.arg(fastk.memory)
  fastk.start <- match.arg(fastk.start)
  qreml.phi.method <- match.arg(qreml.phi.method)

  if (tune.method != "one-step") {
    stop("The optimized engine currently supports tune.method='one-step' only.")
  }

  # dx is constructed inside fgee(), so the optimized engine can update it
  # by reference without making a second full long-data copy. Preserve the
  # original pffr/model-row order before canonicalizing cluster operations;
  # the embedded mgcv model must receive eta/mu in its original row order.
  dx <- data.table::as.data.table(dx)
  original_row_col <- ".fgee_original_row"
  if (original_row_col %in% names(dx)) {
    stop("Reserved column already present in working data: ", original_row_col)
  }
  dx[, (original_row_col) := .I]
  data.table::setorderv(dx, c("cname_", index_long, index_fn))

  fam_key <- .fgee_family_key(glmfit$family)
  link <- tolower(glmfit$family$link)
  use_exact <- isTRUE(exact) && fam_key %in% c("gaussian", "normal") && link == "identity"
  beta0 <- as.numeric(glmfit$coefficients)
  beta_current <- beta0

  if (!all(X_ %in% names(dx))) stop("Some model-matrix columns are missing from dx.")
  if (length(beta0) != length(X_)) stop("Coefficient and model-matrix dimensions differ.")

  method <- .fgee_resolve_sp_method(sp.method, glmfit$family, glmfit$family$link)
  if (method == "legacy") stop("sp.method='legacy' must use working.engine='legacy'.")

  retain <- .fgee_resolve_working_retain(
    requested = working.retain,
    joint.CI = joint.CI,
    var.type = var.type,
    exact = use_exact,
    max.iter = max.iter,
    sp.method = method,
    tune.method = tune.method,
    gee.fit = gee.fit
  )

  fold_info <- NULL
  needs_fastk <- isTRUE(gee.fit) && method %in% c("fastk_staged", "fastk_grad", "qreml_fastk")
  if (needs_fastk) {
    fold_info <- .fgee_make_cluster_folds(as.character(clusters), K = fastk.K,
                                          seed = fastk.seed)
  }

  initial <- .fgee_update_mean_and_corr(
    dx = dx,
    beta = beta0,
    namesd = X_,
    glmfit = glmfit,
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    rho.smooth = rho.smooth,
    rho.pool = rho.pool,
    linpred_method = linpred_method,
    clip_mu = clip_mu,
    update_nuisance = "fixed"
  )
  dx <- initial$dx

  # The value frozen here builds the working variance for every candidate
  # smoothing parameter, so record where it came from, and refuse to continue
  # for the two families whose variance is meaningless without it. This guards
  # against get_family_info()'s silent fallback to 1.0, which for beta yields
  # the entirely plausible-looking and entirely wrong v = mu(1-mu)/2.
  nuisance_frozen <- initial$nuisance
  .fgee_assert_nuisance_frozen(glmfit, nuisance_frozen)

  working0 <- fgee_build_working_stats(
    dx = dx,
    namesd = X_,
    cname_ = "cname_",
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    fpca_fn = initial$fpca,
    retain = retain,
    corr_solver = corr.solver,
    beta0 = beta0,
    exact_gaussian = use_exact,
    gaussian_fold_id = if (fam_key %in% c("gaussian", "normal") && needs_fastk) {
      fold_info$fold_id
    } else NULL,
    copy_dt = FALSE,
    ensure_order = FALSE
  )

  ps <- penalty_setup(glmfit, unpenalized = glmfit$nsdf)

  if (isFALSE(gee.fit)) {
    lambda <- .fgee_initial_lambda(glmfit, ps)
    penalty <- penalty_from_setup(ps, lambda)
    vb <- fgee_var_from_stats(
      working = working0,
      beta = beta0,
      penalty_diag = penalty,
      working0 = working0,
      beta0 = beta0,
      B = boot.samps,
      var.type = var.type,
      exact = FALSE,
      boot.base = "initial"
    )
    glmfit$Vp <- vb
    ci <- compute_joint_ci_compact(
      joint.CI = joint.CI,
      glmfit = glmfit,
      working0 = working0,
      working2 = working0,
      beta0 = beta0,
      beta2 = beta0,
      penalty_diag = penalty,
      bs = bs,
      index = index,
      exact = FALSE
    )
    glmfit <- ci$glmfit
    model_order <- order(dx[[original_row_col]])
    model_eta <- as.numeric(dx$eta[model_order])
    model_fitted <- as.numeric(dx$p[model_order])
    dx[, (original_row_col) := NULL]
    return(list(
      beta = beta0,
      vb = vb,
      rho = NULL,
      di0 = .fgee_compat_scores(working0),
      wi0 = .fgee_compat_bread(working0),
      di = NULL,
      wi = NULL,
      Wbar0 = working0$W_bar,
      Wbar = working0$W_bar,
      working0 = working0,
      working = working0,
      model = glmfit,
      pen.mat = penalty,
      lambda = lambda,
      qn = ci$qn,
      crit = glmfit$crit,
      ci_newdata = glmfit$ci_newdata,
      n_iter = 0L,
      converged = TRUE,
      exact = FALSE,
      sp.method = "initial_frem",
      tuning = NULL,
      working.retain = retain,
      working.engine = "optimized",
      model_eta = model_eta,
      model_fitted = model_fitted,
      data = dx
    ))
  }

  tuning <- fgee_tune_smoothing_optimized(
    working = working0,
    dx = dx,
    namesd = namesd,
    fit.initial = glmfit,
    cv.grid = cv.grid,
    sp.method = method,
    K = if (is.null(fold_info)) fastk.K else fold_info$K,
    seed = fastk.seed,
    sets = if (is.null(fold_info)) NULL else fold_info$sets,
    exact = use_exact,
    fastk.memory = fastk.memory,
    fastk.start = fastk.start,
    qreml.phi.method = qreml.phi.method,
    qreml.phi.fixed = qreml.phi.fixed,
    qreml.phi.weight = qreml.phi.weight,
    qreml.phi.clip = qreml.phi.clip,
    keep.workspace = keep.tuning.workspace,
    verbose = verbose.tuning,
    cname_ = "cname_"
  )

  penalty <- tuning$penalty_mat
  lambda <- tuning$lambda
  coef.diff <- Inf
  gee.itr <- 0L
  working_iter <- working0
  nuisance <- initial$nuisance
  rho_iter <- initial$rho
  fpca_iter <- initial$fpca
  first <- TRUE

  while (coef.diff > beta.tol && gee.itr < max.iter) {
    gee.itr <- gee.itr + 1L

    if (!first) {
      current <- .fgee_update_mean_and_corr(
        dx = dx,
        beta = beta_current,
        namesd = X_,
        glmfit = glmfit,
        corr_fn = corr_fn,
        corr_long = corr_long,
        index_fn = index_fn,
        index_long = index_long,
        rho.smooth = rho.smooth,
        rho.pool = rho.pool,
        linpred_method = linpred_method,
        clip_mu = clip_mu,
        update_nuisance = "moment",
        nuisance = nuisance
      )
      dx <- current$dx
      nuisance <- current$nuisance
      rho_iter <- current$rho
      fpca_iter <- current$fpca
      working_iter <- fgee_build_working_stats(
        dx = dx,
        namesd = X_,
        cname_ = "cname_",
        corr_fn = corr_fn,
        corr_long = corr_long,
        index_fn = index_fn,
        index_long = index_long,
        fpca_fn = fpca_iter,
        retain = retain,
        corr_solver = corr.solver,
        beta0 = beta_current,
        exact_gaussian = use_exact,
        copy_dt = FALSE,
        ensure_order = FALSE
      )
    }
    first <- FALSE

    beta_new <- .fgee_one_step_update(
      working_iter,
      beta = beta_current,
      penalty = penalty,
      exact = use_exact
    )
    coef.diff <- max(abs(beta_new - beta_current))
    beta_current <- beta_new
    glmfit$coefficients <- beta_current
  }

  if (gee.itr >= max.iter && coef.diff > beta.tol && max.iter > 1L) {
    warning("GEE did not converge within max.iter=", max.iter,
            " (final max|dBeta|=", signif(coef.diff, 6), ").")
  }

  beta2 <- as.numeric(beta_current)
  glmfit$coefficients <- beta2
  attr(dx, "nuisance") <- nuisance

  final <- .fgee_update_mean_and_corr(
    dx = dx,
    beta = beta2,
    namesd = X_,
    glmfit = glmfit,
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    rho.smooth = rho.smooth,
    rho.pool = rho.pool,
    linpred_method = linpred_method,
    clip_mu = clip_mu,
    update_nuisance = "moment",
    nuisance = nuisance
  )
  dx <- final$dx

  working2 <- fgee_build_working_stats(
    dx = dx,
    namesd = X_,
    cname_ = "cname_",
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    fpca_fn = final$fpca,
    retain = retain,
    corr_solver = corr.solver,
    beta0 = beta2,
    exact_gaussian = use_exact,
    copy_dt = FALSE,
    ensure_order = FALSE
  )

  vb <- fgee_var_from_stats(
    working = working2,
    beta = beta2,
    penalty_diag = penalty,
    working0 = working0,
    beta0 = beta0,
    B = boot.samps,
    var.type = var.type,
    exact = use_exact,
    boot.base = "initial"
  )
  glmfit$Vp <- vb

  ci <- compute_joint_ci_compact(
    joint.CI = joint.CI,
    glmfit = glmfit,
    working0 = working0,
    working2 = working2,
    beta0 = beta0,
    beta2 = beta2,
    penalty_diag = penalty,
    bs = bs,
    index = index,
    exact = use_exact
  )
  glmfit <- ci$glmfit

  model_order <- order(dx[[original_row_col]])
  model_eta <- as.numeric(dx$eta[model_order])
  model_fitted <- as.numeric(dx$p[model_order])
  dx[, (original_row_col) := NULL]

  list(
    beta = beta2,
    vb = vb,
    rho = final$rho,
    di0 = .fgee_compat_scores(working0),
    wi0 = .fgee_compat_bread(working0),
    di = .fgee_compat_scores(working2),
    wi = .fgee_compat_bread(working2),
    Wbar0 = working0$W_bar,
    Wbar = working2$W_bar,
    working0 = working0,
    working = working2,
    model = glmfit,
    pen.mat = penalty,
    lambda = lambda,
    qn = ci$qn,
    crit = glmfit$crit,
    ci_newdata = glmfit$ci_newdata,
    n_iter = gee.itr,
    converged = (max.iter == 1L) || (coef.diff <= beta.tol),
    exact = use_exact,
    sp.method = method,
    tuning = tuning,
    working.retain = retain,
    working.engine = "optimized",
    model_eta = model_eta,
    model_fitted = model_fitted,
    data = dx
  )
}
