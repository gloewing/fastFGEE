#' Fit a one-step functional generalized estimating equation
#'
#' Fits longitudinal functional-response generalized estimating equations using
#' the validated one-step estimator.  Experimental exact-GLS, pffr-only, legacy
#' engine, and fully iterated controls are intentionally not part of the public
#' interface; their implementations remain internal for regression testing and
#' future method development.
#'
#' @param formula A model formula. The left-hand side should be a functional
#'   response stored as a matrix-like column, typically wrapped in `I()`.
#' @param data A data frame containing the variables in `formula`.
#' @param cluster Name of the cluster identifier column.
#' @param family A family object or family name understood by `refund::pffr()`.
#' @param corr_fn Working correlation in the functional direction: one of
#'   `"independent"`, `"exchangeable"`, `"ar1"`, or `"fpca"`.
#' @param corr_long Working correlation in the longitudinal direction: one of
#'   `"independent"`, `"exchangeable"`, or `"ar1"`. Note that setting exactly
#'   one of `corr_fn`/`corr_long` to `"independent"` forfeits the separable
#'   (Kronecker) structure, so a "simpler" working correlation is not
#'   necessarily cheaper to fit; see `rho.pool`.
#' @param time Optional name of the longitudinal ordering variable.
#' @param long.dir Logical retained for backwards compatibility.
#' @param var.type Variance estimator: `"sandwich"`, `"fastboot"`, or `"boot"`.
#'   `"fastboot"` is preferred to `"boot"` because it is more thoroughly tested.
#' @param pffr.mod Optional fitted `refund::pffr()` object used as the initial
#'   estimator.
#' @param knots Number of spline knots for the initial `pffr()` fit.
#' @param bs Basis type passed to `refund::pffr()`.
#' @param cv Cross-validation mode used for smoothing-parameter selection.
#'   The public one-step interface currently supports only `"fastkfold"`.
#' @param cv.grid Optional grid or staged grids of smoothing parameters.
#' @param rho.smooth Logical; smooth pointwise working-correlation estimates.
#' @param rho.pool Controls pooling of the working-correlation parameter when
#'   only one direction is correlated. `"fn"` (default) estimates a single
#'   pooled functional correlation when `corr_long = "independent"`; `"long"`
#'   pools the longitudinal correlation when `corr_fn = "independent"`;
#'   `"both"` pools either; `"none"` lets the surviving parameter vary over the
#'   opposite index, which reproduces the behaviour of earlier versions. When
#'   both directions are correlated the working correlation is separable and
#'   both parameters are pooled regardless of this setting. Pooling the
#'   functional direction restores that separability, so the correlation
#'   inverse is applied once per cluster rather than once per longitudinal
#'   observation; it is the default because a functional correlation indexed by
#'   visit has no partial-pooling model behind it and is poorly determined when
#'   cluster sizes are unbalanced.
#' @param joint.CI Controls interval construction. `"wild"` requests the
#'   studentized wild-cluster procedure and its optional simultaneous band.
#'   A requested band is a simultaneous-inference procedure, not a guarantee of
#'   nominal finite-sample coverage; reliability depends on the number and
#'   leverage distribution of independent clusters and on nuisance/correlation
#'   stability.
#' @param linpred_method Method used to form linear predictors.
#' @param clip_mu Lower bound used for numerical stabilization of fitted means.
#' @param m.pffr Penalty-order specification passed to `refund::pffr()`.
#' @param check_alignment Logical; check alignment between the wide data and the
#'   long representation produced by `pffr()`.
#' @param boot.samps Number of bootstrap replicates when applicable.
#' @param sp.method Smoothing-parameter selector. `"auto"` uses staged fastK for
#'   Gaussian identity-link models and `"fastk_grad_fast"` otherwise.  The
#'   available explicit selectors are `"fastk_staged"`, `"fastk_grad"`,
#'   `"fastk_grad_fast"`, `"sandwich_qreml"`, and `"qreml_fastk"`.
#'   `"qreml_fastk"` is retained for reproducibility but normally reaches the
#'   same exact fastK solution by a slower route.
#' @param working.retain Retention profile for compact working statistics:
#'   `"auto"`, `"scores"`, `"aggregate"`, or `"full"`.
#' @param corr.solver Correlation inverse backend. `"auto"` uses the package's
#'   direct regular-grid and irregular-grid operators when available;
#'   `"supergauss"` forces the optional Toeplitz backend.
#' @param fastk.K Number of cluster folds for optimized fastK tuning.
#' @param fastk.seed Seed used to construct optimized fastK folds.
#' @param fastk.memory Non-Gaussian fastK workspace: `"balanced"`, `"speed"`,
#'   or deprecated `"lowmem"`.  `"balanced"` is the default and the only layout
#'   eligible for the compiled loss/gradient kernel.  `"lowmem"` is retained
#'   temporarily for backwards compatibility but has not reduced measured peak
#'   process memory and is slated for removal.
#' @param fastk.start Initialization strategy for analytic-gradient fastK.
#' @param fastk.kernel Logical; allow the compiled loss/gradient kernel when the
#'   selected family and the `"balanced"` workspace support it.
#' @param qreml.phi.method Information-scaling diagnostic used by sandwich qREML.
#' @param qreml.phi.fixed Fixed information scale when
#'   `qreml.phi.method = "fixed"`.
#' @param qreml.phi.weight Weight applied to the information mismatch before
#'   clipping.
#' @param qreml.phi.clip Lower and upper safeguards for the qREML information
#'   scale.
#' @param keep.tuning.workspace Logical; retain large tuning workspaces for
#'   debugging.
#' @param verbose.tuning Logical; print optimized tuning progress.
#' @param keep.data Logical; retain long-format working data in the fitted object.
#' @param keep.initial.fit Logical; retain the initial `pffr()` fit.
#' @param keep.working.stats Logical; retain compact initial and final working
#'   statistics.
#' @param ... Reserved only to detect removed estimator controls. Other entries
#'   are rejected rather than silently ignored.
#'
#' @return An object of class `"fgee1step"`. Important components include
#'   `beta`, `vb`, `model`, `lambda`, `tuning`, and, when retained, `working0`
#'   and `working`.
#'
#' @details
#' Public calls always use one coefficient update and the optimized working
#' engine. FastK selectors use the one-step fast cluster cross-validation
#' construction; `"sandwich_qreml"` remains an explicit opt-in selector rather
#' than a cross-validation method. The historical exact, legacy, pffr-only, and
#' fully iterated paths are unexported development interfaces.
#'
#' Pointwise intervals and simultaneous bands answer different inferential
#' questions. A simultaneous band controls a whole-curve maximum statistic in
#' the asymptotic/resampling procedure; it should not be described as guaranteed
#' finite-sample coverage, particularly with very few independent clusters.
#'
#' @references
#' Loewinger, G., Levis, A. W., Cui, E., and Pereira, F. (2025). Fast Penalized
#' Generalized Estimating Equations for Large Longitudinal Functional Datasets.
#'
#' @examples
#' \dontrun{
#' data("d", package = "fastFGEE")
#' fit <- fgee(
#'   formula = Y ~ X1 + X2,
#'   data = d,
#'   cluster = "ID",
#'   family = binomial(link = "logit"),
#'   time = "time",
#'   corr_long = "exchangeable",
#'   corr_fn = "independent"
#' )
#' fgee.plot(fit)
#' }
#'
#' @export
fgee <- function(
    formula,
    data,
    cluster,
    family,
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
    rho.smooth = FALSE,
    rho.pool = c("fn", "none", "long", "both"),
    joint.CI = "wild",
    linpred_method = c("accumulate", "matrix"),
    clip_mu = 0,
    m.pffr = c(2, 1),
    check_alignment = TRUE,
    boot.samps = 3000,
    sp.method = c(
      "auto", "fastk_staged", "fastk_grad", "fastk_grad_fast",
      "sandwich_qreml", "qreml_fastk"
    ),
    working.retain = c("auto", "scores", "aggregate", "full"),
    corr.solver = c("auto", "exact", "supergauss"),
    fastk.K = 10L,
    fastk.seed = 1L,
    fastk.memory = c("balanced", "speed", "lowmem"),
    fastk.start = c("qreml", "robust", "compact"),
    fastk.kernel = fgee_fastk_kernel_ok(),
    qreml.phi.method = c("penalized", "all", "fixed"),
    qreml.phi.fixed = NULL,
    qreml.phi.weight = 1,
    qreml.phi.clip = c(1, 8),
    keep.tuning.workspace = FALSE,
    verbose.tuning = TRUE,
    keep.data = TRUE,
    keep.initial.fit = TRUE,
    keep.working.stats = TRUE,
    ...) {

  public_call <- match.call()
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("", length(dots))
    disabled <- c(
      "exact", "gls", "gee.fit", "gee_fit", "max.iter", "max_iter",
      "niter", "tune.method", "tune_method", "working.engine",
      "working_engine", "fully.iterated", "fully_iterated",
      "full.iterated", "full_iterated", "iterated"
    )
    hit <- which(tolower(dot_names) %in% disabled)
    if (length(hit)) {
      shown <- dot_names[hit]
      shown[!nzchar(shown)] <- "<unnamed>"
      stop(
        "Only the validated one-step estimator is available through fgee(). ",
        "Removed estimator control(s): ", paste(unique(shown), collapse = ", "),
        ". Experimental exact, legacy, and fully iterated engines remain ",
        "internal for package development.",
        call. = FALSE
      )
    }
    shown <- dot_names
    shown[!nzchar(shown)] <- "<unnamed>"
    stop(
      "Unused argument(s) in fgee(): ", paste(unique(shown), collapse = ", "),
      ". The public interface does not silently ignore entries in `...`.",
      call. = FALSE
    )
  }

  if (!is.character(cv) || length(cv) != 1L ||
      !identical(cv, "fastkfold")) {
    stop(
      "The public one-step interface currently supports only cv = ",
      "\"fastkfold\". Other cross-validation and legacy paths remain ",
      "internal for method development.",
      call. = FALSE
    )
  }

  linpred_method <- match.arg(linpred_method)
  sp.method <- match.arg(sp.method)
  working.retain <- match.arg(working.retain)
  corr.solver <- match.arg(corr.solver)
  fastk.memory <- match.arg(fastk.memory)
  fastk.start <- match.arg(fastk.start)
  qreml.phi.method <- match.arg(qreml.phi.method)

  if (identical(fastk.memory, "lowmem")) {
    .Deprecated(
      new = 'fastk.memory = "balanced"',
      package = "fastFGEE",
      msg = paste(
        "fastk.memory = \"lowmem\" is deprecated: validation found no regime",
        "in which it reduced peak process memory, and it cannot use the compiled",
        "fastK kernel. Use \"balanced\" (recommended) or \"speed\"."
      )
    )
  }

  out <- .fgee_fit_internal(
    formula = formula,
    data = data,
    cluster = cluster,
    family = family,
    corr_fn = corr_fn,
    corr_long = corr_long,
    time = time,
    long.dir = long.dir,
    var.type = var.type,
    pffr.mod = pffr.mod,
    knots = knots,
    bs = bs,
    cv = cv,
    cv.grid = cv.grid,
    exact = FALSE,
    rho.smooth = rho.smooth,
    rho.pool = rho.pool,
    joint.CI = joint.CI,
    gee.fit = TRUE,
    linpred_method = linpred_method,
    clip_mu = clip_mu,
    m.pffr = m.pffr,
    check_alignment = check_alignment,
    max.iter = 1L,
    tune.method = "one-step",
    boot.samps = boot.samps,
    working.engine = "optimized",
    sp.method = sp.method,
    working.retain = working.retain,
    corr.solver = corr.solver,
    fastk.K = fastk.K,
    fastk.seed = fastk.seed,
    fastk.memory = fastk.memory,
    fastk.start = fastk.start,
    fastk.kernel = fastk.kernel,
    qreml.phi.method = qreml.phi.method,
    qreml.phi.fixed = qreml.phi.fixed,
    qreml.phi.weight = qreml.phi.weight,
    qreml.phi.clip = qreml.phi.clip,
    keep.tuning.workspace = keep.tuning.workspace,
    verbose.tuning = verbose.tuning,
    keep.data = keep.data,
    keep.initial.fit = keep.initial.fit,
    keep.working.stats = keep.working.stats
  )

  if (is.list(out)) out$call <- public_call
  out
}
