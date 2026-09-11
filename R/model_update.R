#' @keywords internal
#' @noRd
fgee_model_update <- function(mod.fit,
                              MM = NULL,
                              eta = NULL,
                              fitted = NULL,
                              coef_name = "beta",
                              store_in = c("model", "self"),
                              lp_name = "linear.predictors",
                              fitted_name = "fitted.values") {

  store_in <- match.arg(store_in)

  if (is.null(mod.fit[[coef_name]])) {
    stop("mod.fit[['", coef_name, "']] not found (expected coefficient vector).")
  }
  beta <- as.numeric(mod.fit[[coef_name]])

  # The optimized engine already has the final linear predictor and mean in its
  # long working table. Reusing them avoids rebuilding/copying a potentially
  # very large M x p model matrix solely to update the fitted mgcv object.
  if (is.null(eta)) {
    if (is.null(MM)) stop("Provide either MM or eta.")
    eta <- drop(MM %*% beta)
  } else {
    eta <- as.numeric(eta)
  }
  if (!length(eta) || any(!is.finite(eta))) {
    stop("eta must be a non-empty finite vector.")
  }

  tgt <- if (store_in == "model") mod.fit$model else mod.fit
  if (is.null(tgt)) stop("Target model object is NULL; set store_in correctly.")
  if (is.function(tgt)) stop("Target model object is a function, not a fit.")

  tgt$coefficients <- beta
  if (!is.null(tgt$coef) || "coef" %in% names(tgt)) tgt$coef <- beta
  tgt[[lp_name]] <- eta

  if (is.null(fitted)) {
    if (is.null(tgt$family)) {
      stop("Could not determine family. Expected $family on the fitted object.")
    }
    f_link <- gee_family_fns(
      family = tgt$family,
      dispersion = .fgee_nuisance_value(mod.fit, "dispersion", suppressWarnings(mod.fit$rho$dispersion)),
      theta = .fgee_nuisance_value(mod.fit, "theta", suppressWarnings(mod.fit$rho$theta)),
      precision = .fgee_nuisance_value(mod.fit, "precision", suppressWarnings(mod.fit$rho$precision))
    )
    fitted <- f_link$linkinv(eta)
  }
  fitted <- as.numeric(fitted)
  if (length(fitted) != length(eta) || any(!is.finite(fitted))) {
    stop("fitted must be finite and have the same length as eta.")
  }
  tgt[[fitted_name]] <- fitted

  if (store_in == "model") {
    mod.fit$model <- tgt
  } else {
    mod.fit <- tgt
  }

  mod.fit
}
