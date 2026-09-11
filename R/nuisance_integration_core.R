# Integrated nuisance plumbing for fastFGEE -----------------------------------
#
# The underlying implementations have explicit internal names in their defining
# files.  These wrappers add typed nuisance handling without capturing or
# redefining functions according to source-file collation order.

#' @keywords internal
#' @noRd
fgee_update_working_cols_dt <- function(...) {
  args <- .fgee_nuisance_map_dots(.fgee_update_working_cols_dt_core,
                                  list(...))
  family <- args$family
  key <- .fgee_nuisance_family_key(family)
  if (key == "none" || isTRUE(args$exact)) {
    return(do.call(.fgee_update_working_cols_dt_core, args))
  }

  requested <- .fgee_nuisance_or(args$update_nuisance, "fixed")
  requested <- match.arg(requested, c("fixed", "moment", "profile", "auto"))
  param <- .fgee_nuisance_parameter_name(key)
  previous <- args[[param]]
  if (!.fgee_positive_scalar(previous)) {
    from_family <- .fgee_nuisance_from_family(family, key)
    if (!is.null(from_family)) previous <- from_family$value
  }

  # The core updater calculates eta, mu and muprime.  Hold the scalar nuisance
  # fixed for that calculation, then update it once (when requested) and refresh
  # variance and standardized-residual columns from the same typed record.
  core_args <- args
  core_args$update_nuisance <- "fixed"
  if (.fgee_positive_scalar(previous)) {
    core_args[[param]] <- as.numeric(previous)[1L]
  }
  dd <- do.call(.fgee_update_working_cols_dt_core, core_args)

  cluster <- if ("cname_" %in% names(dd)) dd[["cname_"]] else NULL
  control <- getOption("fastFGEE.nuisance.control", list())

  if (requested == "fixed" && .fgee_positive_scalar(previous)) {
    nuisance <- list(
      family = key,
      parameter = param,
      value = as.numeric(previous)[1L],
      method = "fixed",
      source = "supplied_or_model",
      converged = TRUE,
      boundary = FALSE,
      weighting = .fgee_nuisance_or(control$weighting, "observation"),
      n_used = nrow(dd),
      clipped_y = 0L,
      clipped_mu = 0L,
      message = NULL
    )
    nuisance[[param]] <- nuisance$value
    class(nuisance) <- c("fgee_nuisance", "list")
  } else {
    method <- if (requested == "moment") "moment" else "profile"
    if (requested == "fixed" && !.fgee_positive_scalar(previous)) {
      warning(
        "No valid initial ", param, " was supplied for ", key,
        "; estimating it once from the current fitted mean rather than ",
        "silently using 1.",
        call. = FALSE
      )
      method <- "profile"
    }
    nuisance <- fgee_estimate_nuisance(
      y = dd[["Y"]],
      mu = dd[["p"]],
      family = family,
      method = method,
      cluster = cluster,
      start = previous,
      control = control
    )
  }

  .fgee_apply_nuisance_to_working_data(dd, nuisance)
}

#' @keywords internal
#' @noRd
fgee_build_working_stats <- function(...) {
  args <- .fgee_nuisance_map_dots(.fgee_build_working_stats_core,
                                  list(...))
  dx <- args$dx
  out <- do.call(.fgee_build_working_stats_core, args)
  out$nuisance <- attr(dx, "nuisance", exact = TRUE)
  out
}

#' @keywords internal
#' @noRd
get_family_info <- function(glmfit) {
  out <- .fgee_get_family_info_core(glmfit)
  nuisance <- .fgee_nuisance_from_fit(glmfit)
  if (is.null(nuisance) && !is.null(out$nuisance)) nuisance <- out$nuisance
  out$nuisance <- nuisance
  if (!is.null(nuisance)) {
    if (identical(nuisance$family, "beta") ||
        identical(nuisance$family, "negbinomial")) {
      out$dispersion_cpp <- nuisance$value
    }
    if (identical(nuisance$family, "gamma")) {
      out$dispersion <- nuisance$value
    }
  }
  out
}

#' @keywords internal
#' @noRd
.fgee_attach_nuisance_metadata <- function(out) {
  if (!is.list(out)) return(out)

  initial <- if (!is.null(out$working0)) out$working0$nuisance else NULL
  final <- if (!is.null(out$working)) out$working$nuisance else NULL
  if (is.null(final) && !is.null(out$data)) {
    final <- attr(out$data, "nuisance", exact = TRUE)
  }

  existing <- out$nuisance
  if (is.null(existing) || !is.list(existing)) existing <- list()
  if (is.null(existing$initial)) existing$initial <- initial
  if (is.null(existing$final)) existing$final <- final
  existing$tuning <- .fgee_nuisance_or(existing$tuning, "fixed")
  existing$note <- .fgee_nuisance_or(
    existing$note,
    paste(
      "Scalar nuisance parameters are held fixed during smoothing-parameter",
      "tuning and updated only at the designated final working-state update."
    )
  )
  out$nuisance <- existing
  out
}

# Full historical/development engine.  It is deliberately unexported so exact
# Gaussian and fully iterated experiments remain available to package tests and
# method development without being presented as supported user estimators.
#' @keywords internal
#' @noRd
.fgee_fit_internal <- function(...) {
  .fgee_attach_nuisance_metadata(.fgee_fit_internal_core(...))
}
