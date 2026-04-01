#' Fits functional generalized estimating equations (fGEE) for
#' longitudinal functional outcomes using a one-step estimator,
#' with optional fully iterated final estimation.
#'
#' @param formula A model formula. The left-hand side should be a functional
#'   response stored as a matrix-like column (typically wrapped in `I()`).
#' @param data A data frame containing the variables in `formula`.
#' @param cluster Name of the cluster identifier column.
#' @param family A family object or family name understood by `refund::pffr()`.
#' @param corr_fn Working correlation in the functional direction.
#'   One of `"independent"`, `"exchangeable"`, `"ar1"`, or `"fpca"`.
#' @param corr_long Working correlation in the longitudinal direction.
#'   One of `"independent"`, `"exchangeable"`, or `"ar1"`.
#' @param time Optional name of the longitudinal ordering variable.
#' @param long.dir Logical; retained for backwards compatibility.
#' @param var.type Variance estimator. One of `"sandwich"`, `"fastboot"`,
#'   or `"boot"`. For bootstraping, we recommend "fastboot" over "boot" because it is more thoroughly tested.
#' @param pffr.mod Optional fitted `refund::pffr()` object to use as the
#'   initial estimator.
#' @param knots Number of spline knots for the initial `pffr()` fit.
#' @param bs Basis type passed to `refund::pffr()`.
#' @param cv Cross-validation mode used for smoothing parameter selection.
#' @param cv.grid Optional grid (or staged grids) of smoothing parameters.
#' @param exact Logical; for Gaussian identity-link models, use the exact
#'   penalized weighted least-squares update.
#' @param rho.smooth Logical; smooth pointwise correlation estimates over the
#'   relevant index when applicable.
#' @param joint.CI Logical or character controlling joint confidence intervals.
#' @param gee.fit Logical; if `FALSE`, return uncertainty quantification around
#'   the initial fit without the GEE update.
#' @param linpred_method Method used to form linear predictors.
#' @param clip_mu Lower bound used for numerical stabilization of fitted means.
#' @param m.pffr Penalty order specification passed to `refund::pffr()`.
#' @param check_alignment Logical; check alignment between the wide data and the
#'   `pffr()` long representation.
#' @param max.iter Maximum number of GEE iterations. `1` gives the one-step fit.
#' @param tune.method Smoothing-parameter tuning method.
#' @param boot.samps Number of bootstrap replicates used when applicable.
#' @param ... Additional arguments reserved for future use.
#' @return A "fgee1step" object. \code{pffr_initial.fit} contains the initial fit refund::pffr object.
#' \code{vb} contains the variance/covariance matrix for the coefficient estimates (sandwich or bootstrap-based).
#' \code{qn} contains the joint CI quantiles. \code{di} and \
#' \code{wi} are lists of length N, with the updated cluster-specific estimating equation and hessian terms (without the penalty).
#' @author Gabriel Loewinger \email{gloewinger@@gmail.com}
#' @examples
#' \donttest{
#' library(fastFGEE)
#' data("d", package = "fastFGEE")
#' fit <- fgee(
#'  formula = Y ~ X1 + X2,
#'  data = d,
#'  cluster = "ID",
#'  family = binomial(link = "logit"),
#'  time = "time",
#'  corr_long = "exchangeable",
#'  corr_fn = "independent")
#'
#'  fgee.plot(fit)
#'  }
#'
#' @references Gabriel Loewinger, Alex W. Levis, Erjia Cui, and Francisco Pereira. (2025).
#' Fast Penalized Generalized Estimating Equations for Large Longitudinal Functional Datasets. \emph{arXiv:2506.20437}.
#'
#' @export

fgee <- function(formula, data, cluster, family,
                 corr_fn = "ar1",
                 corr_long = "ar1",
                 time = NULL,
                 long.dir = TRUE,
                 var.type = "sandwich",
                 pffr.mod = NULL,
                 knots = NULL,
                 bs = "bs",
                 cv = "fastkfold",
                 cv.grid = NULL,
                 exact = FALSE,
                 rho.smooth = FALSE,
                 joint.CI = "wild",
                 gee.fit = TRUE,
                 linpred_method = c("accumulate", "matrix"),
                 clip_mu = 0,
                 m.pffr = c(2,1),
                 check_alignment = TRUE,
                 max.iter = 1,
                 tune.method = c("one-step", "fully-iterated"),
                 boot.samps = 3000,
                 ...) {

  # === Check for duplicate arguments ===
  call_args <- as.list(match.call())[-1]
  if (any(duplicated(names(call_args)))) {
    dupes <- names(call_args)[duplicated(names(call_args))]
    stop("Duplicate argument(s) in function call: ", paste(unique(dupes), collapse = ", "))
  }

  dots <- list(...)
  if (length(dots) > 0 && !is.null(names(dots))) {
    formal_args <- names(formals(fgee))
    formal_args <- formal_args[formal_args != "..."]
    overlap <- intersect(names(dots), formal_args)
    if (length(overlap) > 0) {
      stop("Argument(s) specified both as named parameter and in ...: ",
           paste(overlap, collapse = ", "))
    }
  }

  linpred_method <- match.arg(linpred_method)

  if (!inherits(formula, "formula")) stop("`formula` must be a formula.")
  if (!is.data.frame(data)) stop("`data` must be a data.frame or data.table.")
  if (!cluster %in% names(data)) stop("cluster variable not found in data: ", cluster)

  if(isFALSE(gee.fit)){
    # sandwich estimator around initial fit
    corr_long <- corr_fn <- "independent"
    exact <- FALSE
  }

  corr_long <- tolower(corr_long)
  corr_fn   <- tolower(corr_fn)
  corr_long <- ifelse(corr_long == "independence", "independent", corr_long)
  corr_fn   <- ifelse(corr_fn   == "independence", "independent", corr_fn)


  if ( !corr_long %in% c("exchangeable", "ar1", "independent")) {
    stop("corr_long must be one of 'exchangeable','ar1','independence'.")
  }
  if ( !corr_fn %in% c("exchangeable", "ar1", "independent", "fpca")) {
    stop("corr_fn must be one of 'exchangeable','ar1','independence', 'fpca'.")
  }

  # Decide which index is used to estimate rho(s)
  index <- if (isTRUE(long.dir)) "yindex.vec" else "time"
  if (!isTRUE(long.dir) && is.null(time)) {
    stop("If long.dir=FALSE, you must provide `time`.")
  }

  # Initial pffr fit
  if (is.null(pffr.mod)) {
    Y_nm <- all.vars(formula)[1]
    L <- ncol(data[[Y_nm]])
    if (is.null(knots)) knots <- round(L / 4)

    # Force evaluation of the family argument before passing it to pffr
    fam <- if (is.character(family)) {
      match.fun(family)() # If user passes "Gamma", call Gamma()
    } else {
      family # If user passes Gamma(link="log"), it's already an object
    }

    fit_pffr <- refund::pffr(
      formula = formula,
      algorithm = "bam",
      family = fam,
      discrete = TRUE,
      bs.yindex = list(bs = bs,
                       k = knots+1,
                       m = m.pffr),
      bs.int = list(bs = bs,
                    k = knots+1,
                    m = m.pffr),
      data = data
    )
  } else {
    fit_pffr <- pffr.mod
  }

  # Build long dt + model matrix columns
  built <- .fgee_build_long_data(
    data = data,
    fit_pffr = fit_pffr,
    formula = formula,
    cluster = cluster,
    time = time,
    check_alignment = check_alignment
  )
  dx <- built$dx
  X_cols <- built$X_

  clusters <- unique(dx$cname_)
  N_clusters <- length(clusters)

  if (isTRUE(exact) & fit_pffr$family$family != "gaussian") {
    message("exact=TRUE only valid for family='gaussian'. Using one-step (exact=FALSE)")
    exact <- FALSE
  }


  if(is.null(cv.grid)){
    cv.grid <- vector(length = 3, "list")
    # pre-grid (never include 0 in pre-grid!)
    cv.grid[[1]] <- sort(unique( c(10^(-6:4),
                                   10^(-6:4) * 0.25,
                                   10^(-6:4) * 0.5,
                                   10^(-6:4) * 0.75,
                                   10^(-6:4) * 0.9)))

    cv.grid[[2]] <- c(1e-2, 0.1, 1, 10, 100)
    cv.grid[[3]] <- c(0.5, 0.75, 1, 1.3, 2, 5)
  }

  mod.fit <- fun.gee1step.dist_itr(
    orig.data = dx,
    dx = dx,
    formula = formula,
    # family = family,
    X_ = X_cols,
    Y_ = all.vars(formula)[1],
    namesd = X_cols,
    N_clusters = N_clusters,
    clusters = clusters,
    glmfit = fit_pffr,
    yindex.vec = dx$yindex.vec,
    time = time,
    index = index,
    bs = bs,
    corr_fn = corr_fn,
    corr_long = corr_long,
    #parallel = parallel,
    var.type = var.type,
    cv = cv,
    cv.grid = cv.grid,
    rho.smooth = rho.smooth,
    gee.fit = gee.fit,
    joint.CI = joint.CI,
    max.iter = max.iter,
    tune.method = tune.method,
    exact = exact,
    # V.inv = V.inv,
    linpred_method = linpred_method,
    clip_mu = clip_mu,
    boot.samps = boot.samps
  )

  # update mgcv model object
  MM <- suppressWarnings(stats::model.matrix(fit_pffr))
  mod.fit <- fgee_model_update(mod.fit = mod.fit, MM = MM)

  result <- append(mod.fit, list(
    call = match.call(),
    formula = formula,
    family = family,
    outcome = all.vars(formula)[1],
    xnames = X_cols,
    data = dx,
    pffr_initial.fit = fit_pffr,
    cluster_sizes = as.vector(dx[, .N, keyby = cname_][, N])
  ))
  class(result) <- "fgee1step"
  result
}



#' Internal GEE engine used by \code{\link{fgee}}
#' @keywords internal
#' @noRd
fun.gee1step.dist_itr <- function(orig.data, dx, formula, X_, Y_, namesd,
                                  N_clusters, clusters, glmfit, yindex.vec, bs,
                                  corr_fn = "independent",
                                  corr_long = "independent",
                                  index_fn = "yindex.vec",
                                  index_long = "time",
                                  parallel = FALSE,
                                  var.type = "sandwich",
                                  cv = "one-step",
                                  cv.grid,
                                  boot.samps = 3000,
                                  rho.smooth = TRUE,
                                  joint.CI = TRUE,
                                  time = NULL,
                                  index = "yindex.vec",
                                  folds.list = NULL,
                                  sets = NULL,
                                  preTn = TRUE,
                                  linpred_method = c("accumulate", "matrix"),
                                  clip_mu = 0,
                                  beta.tol = 1e-6,
                                  max.iter = 1,
                                  tune.method = c("one-step", "fully-iterated"),
                                  tune.full_fn = NULL,
                                  exact = FALSE,
                                  gee.fit = TRUE,
                                  eval_prop = 1) {

  linpred_method <- match.arg(linpred_method)
  tune.method <- match.arg(tune.method)

  # ---- isolate from caller + canonical order once
  dx <- data.table::copy(data.table::as.data.table(dx))
  data.table::setorderv(dx, c("cname_", index_long, index_fn))
  dr <- data.table::copy(dx)  # final recompute base

  # ------------------------------------------------------------
  # Exact Gaussian mode flag (minimal integration)
  #   - only meaningful for Gaussian + identity link
  #   - IMPORTANT: we DO NOT set exact=TRUE inside .fgee_update_mean_parts(),
  #     because corr.est must use residuals (Y - mu)/sqrt(v), not Y/sqrt(v).
  # ------------------------------------------------------------
  fam_lc  <- tolower(glmfit$family$family)
  link_lc <- tolower(glmfit$family$link)
  use_exact_gaussian <- isTRUE(exact) && (fam_lc == "gaussian") && (link_lc == "identity")

  # ------------------------------------------------------------
  # Initial estimate (PIVOT beta0)
  # ------------------------------------------------------------
  beta0 <- as.numeric(glmfit$coefficients)
  beta_current <- beta0

  # sanity
  if (!all(X_ %in% names(dx))) stop("Some X_ columns not found in dx.")
  if (length(beta_current) != length(X_)) {
    stop("length(beta_current) != length(X_). You likely need to pass the full design column set in X_.\n",
         "length(beta_current)=", length(beta_current), " length(X_)=", length(X_))
  }

  # ------------------------------------------------------------
  # Mean/variance/residuals at beta0
  #   NOTE: exact is intentionally FALSE here even if use_exact_gaussian=TRUE
  #         so that resid=(Y - p)/sqrtv for corr.est.
  # ------------------------------------------------------------
  fi <- get_family_info(glmfit)

  dx <- fgee_update_working_cols_dt(
    dx = dx,
    namesd = X_,
    beta = beta_current,
    family = glmfit$family,
    link = glmfit$family$link,
    exact = FALSE,                  # ALWAYS score residuals for corr.est
    gaussian_v_by = index_fn,
    update_nuisance = "fixed",     # use this to start out with mgcv:gam() theta/phi parameter estimates
    theta = if (fi$family_cpp == "negbinomial") fi$dispersion_cpp else NULL,
    precision = if (fi$family_cpp == "beta") fi$dispersion_cpp else NULL,
    clamp_eps = max(1e-8, min(clip_mu, 0.499999)),            # ensure no 0/1 exact issues
    linpred_method = linpred_method,
    eta_by = NULL,
    copy = FALSE
  )
  # ------------------------------------------------------------
  # Correlation (+ FPCA if needed) at beta0
  # ------------------------------------------------------------
  cor0 <- corr.est(
    dx = dx,
    cname_ = "cname_",
    index_fn = index_fn,
    index_long = index_long,
    corr_fn = corr_fn,
    corr_long = corr_long,
    resid_col = "resid",
    rho.smooth = rho.smooth,
    ar = "mom",
    glmfit = if (isTRUE(rho.smooth)) glmfit else NULL,
    fpca_fn = NULL,
    copy_dt = FALSE
  )

  dx <- cor0$dx # cor0$rho$rho_fn$rho_fn <- mean(cor0$rho$rho_fn$rho_fn)
  # ------------------------------------------------------------
  # W/D at beta0 (SCORE-based; used for tuning and wild bootstrap pivot)
  # ------------------------------------------------------------
  wi0 <- W.estimate(
    dx,
    namesd = X_,
    cname_ = "cname_",
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    algo = "gschur",
    fpca_fn = cor0$fpca,
    id.vec = clusters,
    copy_dt = FALSE
  )

  gc()

  di0 <- D.estimate(
    dx,
    namesd = X_,
    cname_ = "cname_",
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    algo = "gschur",
    fpca_fn = cor0$fpca,
    id.vec = clusters,
    copy_dt = FALSE
  )
  gc()
  # ------------------------------------------------------------
  # Exact Gaussian RHS at beta0 (ONLY for coefficient updates)
  #   - does NOT replace di0 (di0 stays score-based)
  #   - resid_exact = Y/sqrtv is used ONLY inside D.estimate via resid_col
  # ------------------------------------------------------------
  di0_exact <- NULL
  if (use_exact_gaussian) {

    if (!("Y" %in% names(dx))) stop("Exact Gaussian mode requires column 'Y' in dx.")

    # sqrtv should be present because fgee_update_working_cols_dt() writes it,
    # but make this robust anyway:
    if (!("sqrtv" %in% names(dx))) {
      if ("v" %in% names(dx)) {
        dx[, sqrtv := sqrt(v)]
      } else {
        stop("Exact Gaussian mode requires column 'sqrtv' or 'v' in dx. ",
             "Call fgee_update_working_cols_dt() before forming di_exact.")
      }
    }

    dx[, resid_exact := Y / sqrtv]

    di0_exact <- D.estimate(
      dx,
      namesd = X_,
      cname_ = "cname_",
      corr_fn = corr_fn,
      corr_long = corr_long,
      index_fn = index_fn,
      index_long = index_long,
      algo = "gschur",
      fpca_fn = cor0$fpca,
      id.vec = clusters,
      copy_dt = FALSE,
      resid_col = "resid_exact"
    )
    gc()
    dx[, resid_exact := NULL]
  }

  if (isFALSE(gee.fit)) {
    # only sandwich estimator

    # penalty from pffr()
    pen_setup <- penalty_setup(glmfit, unpenalized = glmfit$nsdf)
    sp_vec <- glmfit[["sp"]]
    if (is.list(sp_vec)) sp_vec <- sp_vec[[1]]
    sp_vec <- as.numeric(sp_vec)
    if (length(sp_vec) == 0L) sp_vec <- 0
    penalty_diag <- penalty_from_setup(pen_setup, lambda = sp_vec)

    vb <- var.est(di = di0, wi = wi0, beta2 = beta0,
                  wi0 = wi0, di0 = di0, beta0 = beta0,
                  penalty_diag = penalty_diag,
                  B = boot.samps,
                  var.type = var.type,
                  exact = FALSE,
                  boot.base = "initial", #c("initial", "final"),
                  return.boot = FALSE,
                  seed = NULL,
                  verbose = FALSE,
                  block_size = NULL # memory-safe fastboot
    )

    glmfit$Vp <- vb

    # ------------------------------------------------------------
    # joint CIs
    #   - wild: pivot at beta0 using di0/wi0, center at beta2
    # ------------------------------------------------------------
    ci_result <- compute_joint_ci(
      joint.CI = joint.CI,
      glmfit = glmfit,
      di0 = di0,
      wi0 = wi0,
      di2 = di0,
      wi2 = wi0,
      beta0 = beta0,
      beta2 = beta0,
      penalty_diag = penalty_diag,
      bs = bs,
      index = index,
      alpha = 0.05,  # Could make this a parameter
      wild_progress_every = 0,
      exact = use_exact_gaussian
    )

    qn <- ci_result$qn
    glmfit <- ci_result$glmfit

    # ------------------------------------------------------------
    # return (unchanged fields + iteration info)
    # ------------------------------------------------------------
    result <- list(
      beta = as.vector(beta0),
      vb = vb,
      rho = NULL,
      di0 = di0,
      wi0 = wi0,
      di = NULL,
      wi = NULL,
      model = glmfit,
      pen.mat = penalty_diag,
      lambda = sp_vec,
      qn = qn,
      crit = glmfit$crit,
      ci_newdata = glmfit$ci_newdata,
      n_iter = 0,
      converged = NULL,
      exact = use_exact_gaussian
    )

    return(result)
  }

  # ------------------------------------------------------------
  # Penalty setup and tuning
  # ------------------------------------------------------------
  # ------------------------------------------------------------
  # Tuning (default: one-step tuning using wi0/di0)
  #   NOTE: For exact Gaussian mode, this still uses SCORE di0.
  #         You said you will update CV for exact later; this avoids
  #         changing CV behavior now.
  # ------------------------------------------------------------
  tuning_result <- tune_smoothing_parameters(
    glmfit = glmfit,
    cv.grid = cv.grid,
    wi0 = wi0,
    di0 = di0,
    di0_exact = di0_exact,
    use_exact_gaussian = use_exact_gaussian,
    namesd = namesd,
    dx = dx,
    tune.method = tune.method,
    tune.full_fn = tune.full_fn,
    folds.list = folds.list,
    sets = sets,
    preTn = preTn,
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    rho.smooth = rho.smooth,
    X_ = X_,
    cv = cv,
    fid_ = index_long,
    eval_prop = eval_prop
  )

  penalty_diag <- tuning_result$penalty_diag
  cv.lambda <- tuning_result$cv.lambda
  # ------------------------------------------------------------
  # ITERATIVE penalized GEE update
  #   max.iter=1 => one-step special case
  #
  # Exact Gaussian mode:
  #   - corr.est uses score residuals as usual
  #   - beta update uses D_exact = X^T V^{-1} Y (via resid_exact column)
  # ------------------------------------------------------------
  coef.diff <- Inf
  gee.itr <- 0

  # reuse first-iteration quantities
  wi_iter <- wi0
  di_iter <- di0
  di_exact_iter <- if (use_exact_gaussian) di0_exact else NULL
  first_iter <- TRUE

  while (coef.diff > beta.tol && gee.itr < max.iter) {

    gee.itr <- gee.itr + 1
    if (max.iter > 1) message("GEE Iteration ", gee.itr)

    if (!first_iter) {

      dx <- fgee_update_working_cols_dt(
        dx = dx,
        namesd = X_,
        beta = beta_current,
        family = glmfit$family,
        link = glmfit$family$link,
        exact = FALSE,                  # IMPORTANT: score residuals for corr.est
        gaussian_v_by = index_fn,
        update_nuisance = "moment",     # or "fixed"
        theta = attributes(dx)$nuisance$theta,
        precision = attributes(dx)$nuisance$precision,
        dispersion =  attributes(dx)$nuisance$dispersion,
        zi_prob = attributes(dx)$nuisance$zi_prob,
        clamp_eps = max(1e-8, min(clip_mu, 0.499999)),
        linpred_method = linpred_method,
        eta_by = NULL,
        copy = FALSE
      )

      cor_iter <- corr.est(
        dx = dx,
        cname_ = "cname_",
        index_fn = index_fn,
        index_long = index_long,
        corr_fn = corr_fn,
        corr_long = corr_long,
        resid_col = "resid",
        rho.smooth = rho.smooth,
        ar = "mom",
        glmfit = if (isTRUE(rho.smooth)) glmfit else NULL,
        fpca_fn = NULL,
        copy_dt = FALSE
      )

      dx <- cor_iter$dx

      wi_iter <- W.estimate(
        dx,
        namesd = X_,
        cname_ = "cname_",
        corr_fn = corr_fn,
        corr_long = corr_long,
        index_fn = index_fn,
        index_long = index_long,
        algo = "gschur",
        fpca_fn = cor_iter$fpca,
        id.vec = clusters,
        copy_dt = FALSE
      )

      # score-based D (kept for consistency / potential later use)
      di_iter <- D.estimate(
        dx,
        namesd = X_,
        cname_ = "cname_",
        corr_fn = corr_fn,
        corr_long = corr_long,
        index_fn = index_fn,
        index_long = index_long,
        algo = "gschur",
        fpca_fn = cor_iter$fpca,
        id.vec = clusters,
        copy_dt = FALSE
      )

      # exact RHS D (only if enabled)
      if (use_exact_gaussian) {

        if (!("sqrtv" %in% names(dx))) {
          if ("v" %in% names(dx)) {
            dx[, sqrtv := sqrt(v)]
          } else {
            stop("Exact Gaussian mode requires column 'sqrtv' or 'v' in dx. ",
                 "Call fgee_update_working_cols_dt() before forming di_exact.")
          }
        }

        dx[, resid_exact := Y / sqrtv]

        di_exact_iter <- D.estimate(
          dx,
          namesd = X_,
          cname_ = "cname_",
          corr_fn = corr_fn,
          corr_long = corr_long,
          index_fn = index_fn,
          index_long = index_long,
          algo = "gschur",
          fpca_fn = cor_iter$fpca,
          id.vec = clusters,
          copy_dt = FALSE,
          resid_col = "resid_exact"
        )

        dx[, resid_exact := NULL]
      }

    }

    first_iter <- FALSE

    wi_bar <- Reduce("+", wi_iter) / N_clusters

    if (use_exact_gaussian) {

      wi_sum <- Reduce("+", wi_iter)
      di_exact_sum <- Reduce("+", di_exact_iter)

      beta_new <- as.numeric(solve_pd(
        a = wi_sum + penalty_diag,
        b = di_exact_sum
      ))


    } else {

      di_sum <- Reduce("+", di_iter)
      penalty_vec <- as.numeric(penalty_diag %*% beta_current)

      step <- solve_pd(
        a = wi_bar + penalty_diag,
        b = di_sum - penalty_vec * N_clusters
      )

      beta_new <- beta_current + as.numeric(step) / N_clusters
    }

    coef.diff <- max(abs(beta_new - beta_current))
    beta_current <- beta_new
    glmfit$coefficients <- beta_current

    if (max.iter > 1) message("  max|dBeta| = ", signif(coef.diff, 6))
  }

  if (max.iter > 1 && gee.itr >= max.iter && coef.diff > beta.tol) {
    warning("GEE did not converge within max.iter=", max.iter,
            " (final max|dBeta|=", signif(coef.diff, 6), ").")
  }

  beta2 <- as.numeric(beta_current)
  glmfit$coefficients <- beta2

  attributes(dr)$nuisance <- attributes(dx)$nuisance # update nuisance parameters

  # free dx if you want; we use dr for final recompute
  rm(dx)

  # ------------------------------------------------------------
  # Recompute mean/variance/residuals at final beta (score residuals)
  # ------------------------------------------------------------
  dr <- fgee_update_working_cols_dt(
    dx = dr,
    namesd = X_,
    beta = beta2,
    family = glmfit$family,
    link = glmfit$family$link,
    exact = FALSE,                  # core residuals for corr.est (always FALSE for sandiwch)
    gaussian_v_by = index_fn,
    update_nuisance = "moment",     # or "fixed"
    dispersion = attributes(dr)$nuisance$dispersion,
    theta = attributes(dr)$theta,
    precision = attributes(dr)$precision,
    zi_prob = attributes(dr)$zi_prob,
    clamp_eps = max(1e-8, min(clip_mu, 0.499999)),
    linpred_method = linpred_method,
    eta_by = NULL,
    copy = FALSE
  )


  # ------------------------------------------------------------
  # Final correlation (+ FPCA) at beta2 (REQUIRED before wi2/di2)
  # ------------------------------------------------------------
  cor2 <- corr.est(
    dx = dr,
    cname_ = "cname_",
    index_fn = index_fn,
    index_long = index_long,
    corr_fn = corr_fn,
    corr_long = corr_long,
    resid_col = "resid",
    rho.smooth = rho.smooth,
    ar = "mom",
    glmfit = if (isTRUE(rho.smooth)) glmfit else NULL,
    fpca_fn = NULL,
    copy_dt = FALSE
  )
  dr <- cor2$dx

  # ------------------------------------------------------------
  # Final W/D at beta2 (score-based D for sandwich)
  # ------------------------------------------------------------
  wi2 <- W.estimate(
    dr,
    namesd = X_,
    cname_ = "cname_",
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    algo = "gschur",
    fpca_fn = cor2$fpca,
    id.vec = clusters,
    copy_dt = FALSE
  )
  gc()
  # for sandwich variance estimator, do not use exact D even when exact = TRUE
  di2 <- D.estimate(
    dr,
    namesd = X_,
    cname_ = "cname_",
    corr_fn = corr_fn,
    corr_long = corr_long,
    index_fn = index_fn,
    index_long = index_long,
    algo = "gschur",
    fpca_fn = cor2$fpca,
    id.vec = clusters,
    copy_dt = FALSE
  )
  gc()
  # ------------------------------------------------------------
  # Variance estimator (uses di2/wi2)
  # ------------------------------------------------------------
  vb <- var.est(di = di2, wi = wi2, beta2 = beta2,
                wi0 = wi0, di0 = di0, beta0 = beta0,
                penalty_diag = penalty_diag,
                B = boot.samps,
                var.type = var.type,
                exact = exact,
                boot.base = "initial", #c("initial", "final"),
                return.boot = FALSE,
                seed = NULL,
                verbose = FALSE,
                block_size = NULL # memory-safe fastboot
  )

  glmfit$Vp <- vb

  # ------------------------------------------------------------
  # joint CIs
  #   - wild: pivot at beta0 using di0/wi0, center at beta2
  # ------------------------------------------------------------
  ci_result <- compute_joint_ci(
    joint.CI = joint.CI,
    glmfit = glmfit,
    di0 = di0,
    wi0 = wi0,
    di2 = di2,
    wi2 = wi2,
    beta0 = beta0,
    beta2 = beta2,
    penalty_diag = penalty_diag,
    bs = bs,
    index = index,
    alpha = 0.05,  # Could make this a parameter
    wild_progress_every = 0,
    exact = use_exact_gaussian
  )

  qn <- ci_result$qn
  glmfit <- ci_result$glmfit

  # ------------------------------------------------------------
  # return (unchanged fields + iteration info)
  # ------------------------------------------------------------
  result <- list(
    beta = as.vector(beta2),
    vb = vb,
    rho = cor2$rho,
    di0 = di0,
    wi0 = wi0,
    di = di2,
    wi = wi2,
    model = glmfit,
    pen.mat = penalty_diag,
    lambda = if (!is.null(cv.lambda)) cv.lambda$lambda.star else NULL,
    qn = qn,
    crit = glmfit$crit,
    ci_newdata = glmfit$ci_newdata,
    n_iter = gee.itr,
    converged = (max.iter == 1L) || (coef.diff <= beta.tol),
    exact = use_exact_gaussian
  )

  return(result)
}

