# Unified smoothing-parameter dispatch ----------------------------------------

.fgee_resolve_sp_method <- function(sp.method, family, link) {
  choices <- c(
    "auto", "fastk_staged", "fastk_grad", "fastk_grad_fast",
    "sandwich_qreml", "qreml_fastk", "legacy"
  )
  sp.method <- match.arg(sp.method, choices)
  if (sp.method != "auto") return(sp.method)

  key <- .fgee_family_key(family)
  if (key %in% c("gaussian", "normal") && tolower(link) == "identity") {
    "fastk_staged"
  } else {
    "fastk_grad_fast"
  }
}

.fgee_qreml_information_summary <- function(info) {
  if (is.null(info)) return(NULL)
  keep <- c(
    "phi_all", "phi_pen", "phi_component", "eig_pen_q10",
    "eig_pen_median", "eig_pen_q90", "anisotropy", "effective_rank",
    "penalized_rank", "H_ridge", "center_scores"
  )
  info[intersect(keep, names(info))]
}

.fgee_compact_tuning_diagnostics <- function(x) {
  method <- .fgee_or(x$method, "unknown")
  common <- list(
    method = method,
    score = .fgee_or(x$score, NA_real_),
    convergence = .fgee_or(x$convergence, NA_integer_),
    message = .fgee_or(x$message, NULL),
    boundary = isTRUE(x$boundary),
    evaluations = .fgee_or(x$evaluations, NA_integer_),
    probe_evaluations = .fgee_or(x$probe_evaluations, NA_integer_),
    gradient_evaluations = .fgee_or(x$gradient_evaluations, NA_integer_),
    log10_multiplier = .fgee_or(x$log10_multiplier, NULL)
  )

  if (grepl("qreml", method, fixed = TRUE) || !is.null(x$phi_used)) {
    common$phi_method <- .fgee_or(x$phi_method, NULL)
    common$phi_raw <- .fgee_or(x$phi_raw, NA_real_)
    common$phi_used <- .fgee_or(x$phi_used, NA_real_)
    common$working_information_scale <- .fgee_or(x$n_eff, NA_real_)
    common$exact_gaussian <- isTRUE(.fgee_or(x$prep$exact, FALSE))
    common$penalty_scale <- .fgee_or(x$prep$penalty_scale, 1)
    common$information <- .fgee_qreml_information_summary(x$information)
  }

  if (identical(method, "fastk_staged")) {
    common$stage_scores <- .fgee_or(x$mse.ls, NULL)
    common$stage_standard_errors <- .fgee_or(x$se.ls, NULL)
    common$final_grid <- .fgee_or(x$grid, NULL)
    common$argmin <- .fgee_or(x$argmin, NA_integer_)
  }

  common
}

.fgee_tune_result <- function(x, method = x$method, keep.workspace = FALSE) {
  lambda <- as.numeric(.fgee_or(x$lambda, x$lambda.star))
  if (!length(lambda)) stop("Tuning result did not contain lambda values.")
  list(
    method = method,
    lambda = lambda,
    lambda.star = matrix(lambda, nrow = 1L),
    penalty_mat = x$penalty_mat,
    score = .fgee_or(x$score, NA_real_),
    convergence = .fgee_or(x$convergence, NA_integer_),
    boundary = isTRUE(x$boundary),
    evaluations = .fgee_or(x$evaluations, NA_integer_),
    diagnostics = .fgee_compact_tuning_diagnostics(x),
    workspace = if (isTRUE(keep.workspace)) .fgee_or(x$workspace, x$prep) else NULL,
    raw = if (isTRUE(keep.workspace)) x else NULL
  )
}

#' Tune smoothing parameters using the optimized working-statistics engine
#' @keywords internal
#' @noRd
fgee_tune_smoothing_optimized <- function(
    working,
    dx,
    namesd,
    fit.initial,
    cv.grid,
    sp.method = c(
      "auto", "fastk_staged", "fastk_grad", "fastk_grad_fast",
      "sandwich_qreml", "qreml_fastk", "legacy"
    ),
    K = 10L,
    seed = 1L,
    sets = NULL,
    exact = FALSE,
    fastk.memory = c("balanced", "speed", "lowmem"),
    fastk.start = c("qreml", "robust", "compact"),
    qreml.phi.method = c("penalized", "all", "fixed"),
    qreml.phi.fixed = NULL,
    qreml.phi.weight = 1,
    qreml.phi.clip = c(1, 8),
    qreml.bound = 6,
    grad.bound = 8,
    fastk.kernel = fgee_fastk_kernel_ok(),
    keep.workspace = FALSE,
    verbose = TRUE,
    cname_ = "cname_") {

  fastk.memory <- match.arg(fastk.memory)
  fastk.start <- match.arg(fastk.start)
  qreml.phi.method <- match.arg(qreml.phi.method)
  method <- .fgee_resolve_sp_method(
    sp.method,
    family = fit.initial$family,
    link = fit.initial$family$link
  )

  if (method == "legacy") {
    stop("Legacy tuning is handled by the legacy fGEE engine.")
  }

  qprep <- NULL
  qfit <- NULL
  need_qreml <- method %in% c("sandwich_qreml", "qreml_fastk") ||
    (method %in% c("fastk_grad", "fastk_grad_fast") && fastk.start == "qreml")

  if (need_qreml) {
    qfit <- tryCatch({
      qprep <- fgee_qreml_prepare(working, fit.initial, exact = exact)
      fgee_tune_qreml(
        qprep,
        phi_method = qreml.phi.method,
        phi_fixed = qreml.phi.fixed,
        phi_weight = qreml.phi.weight,
        phi_clip = qreml.phi.clip,
        bound = qreml.bound,
        verbose = verbose
      )
    }, error = function(e) {
      if (identical(method, "sandwich_qreml")) stop(e)
      warning(
        "The qREML initializer failed; continuing with exact fastK starts. ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    })
  }

  if (method == "sandwich_qreml") {
    upper_hit <- !is.null(qreml.phi.clip) &&
      isTRUE(all.equal(qfit$phi_used, max(qreml.phi.clip)))
    if (isTRUE(qfit$boundary) || upper_hit) {
      warning(
        "sandwich_qreml reached a smoothing or information-scale safeguard. ",
        "Inspect fit$tuning$diagnostics, or use sp.method='qreml_fastk' ",
        "to finish with the exact fastK criterion.",
        call. = FALSE
      )
    }
    return(.fgee_tune_result(qfit, method = method,
                             keep.workspace = keep.workspace))
  }

  fprep <- fgee_fastk_prepare(
    working = working,
    dx = dx,
    namesd = namesd,
    fit.initial = fit.initial,
    K = K,
    seed = seed,
    sets = sets,
    exact = exact,
    memory = fastk.memory,
    cname_ = cname_
  )

  if (method == "fastk_staged") {
    out <- fgee_tune_fastk_staged(fprep, cv.grid = cv.grid, verbose = verbose)
    return(.fgee_tune_result(out, method = method,
                             keep.workspace = keep.workspace))
  }

  if (method == "fastk_grad_fast") {
    out <- fgee_tune_fastk_grad_fast(
      fprep,
      qreml_lambda = if (!is.null(qfit)) qfit$lambda else NULL,
      bound = grad.bound,
      kernel = fastk.kernel,
      verbose = verbose
    )
    ans <- .fgee_tune_result(out, method = method,
                             keep.workspace = keep.workspace)
    if (!is.null(qfit)) {
      ans$diagnostics$qreml_start <- .fgee_compact_tuning_diagnostics(qfit)
    }
    return(ans)
  }

  if (method == "fastk_grad") {
    out <- fgee_tune_fastk_grad(
      fprep,
      qreml_lambda = if (!is.null(qfit)) qfit$lambda else NULL,
      start_strategy = fastk.start,
      bound = grad.bound,
      verbose = verbose
    )
    ans <- .fgee_tune_result(out, method = method,
                             keep.workspace = keep.workspace)
    if (!is.null(qfit)) {
      ans$diagnostics$qreml_start <- .fgee_compact_tuning_diagnostics(qfit)
    }
    return(ans)
  }

  if (method == "qreml_fastk") {
    # The returned lambda always minimizes the exact fastK criterion from the
    # tested starts.  qREML is only a deterministic coefficient-specific start,
    # so working-correlation misspecification can affect speed but not the final
    # criterion being optimized.
    out <- fgee_tune_fastk_grad(
      fprep,
      qreml_lambda = if (is.null(qfit)) NULL else qfit$lambda,
      start_strategy = "qreml",
      bound = grad.bound,
      verbose = verbose
    )
    ans <- .fgee_tune_result(out, method = method,
                             keep.workspace = keep.workspace)
    ans$diagnostics$qreml <- if (is.null(qfit)) {
      list(method = "sandwich_qreml", failed = TRUE)
    } else {
      .fgee_compact_tuning_diagnostics(qfit)
    }
    ans$diagnostics$fastk <- .fgee_compact_tuning_diagnostics(out)
    upper_hit <- !is.null(qfit) && !is.null(qreml.phi.clip) &&
      isTRUE(all.equal(qfit$phi_used, max(qreml.phi.clip)))
    ans$diagnostics$recommend_fastk <- is.null(qfit) ||
      isTRUE(qfit$boundary) || upper_hit
    if (isTRUE(keep.workspace)) {
      ans$raw <- list(qreml = qfit, fastk = out)
      ans$workspace <- list(qreml = qprep, fastk = fprep)
    }
    return(ans)
  }

  stop("Unhandled sp.method='", method, "'.")
}
