#' @keywords internal
#' @noRd
fgee_model_update <- function(mod.fit, MM,
                              coef_name = "beta",
                              store_in = c("model", "self"),
                              lp_name = "linear.predictors",
                              fitted_name = "fitted.values") {

  store_in <- match.arg(store_in)

  if (is.null(mod.fit[[coef_name]])) {
    stop("mod.fit[['", coef_name, "']] not found (expected coefficient vector).")
  }
  beta <- as.numeric(mod.fit[[coef_name]])

  if (missing(MM) || is.null(MM)) stop("MM must be provided (design matrix).")

  # allow Matrix objects etc.
  eta <- drop(MM %*% beta)

  # Where is the fitted model object stored?
  tgt <- if (store_in == "model") mod.fit$model else mod.fit

  if (is.null(tgt)) stop("Target model object is NULL; set store_in correctly.")
  if (is.function(tgt)) stop("Target model object is a function, not a fit.")

  # ---- set coefficients
  tgt$coefficients <- beta
  if (!is.null(tgt$coef) || "coef" %in% names(tgt)) tgt$coef <- beta

  # ---- linear predictor
  tgt[[lp_name]] <- eta

  # robust link/family extraction is already done by gee_family_fns,
  # we just need to pass it the family object.
  if (is.null(tgt$family)) stop("Could not determine family. Expected $family on the fitted object.")

  # 1. Call the modern helper function to get the correct linkinv function
  f_link <- gee_family_fns(
    family = tgt$family,
    # Pass nuisance parameters from the final model fit if they exist.
    # Use suppressWarnings to avoid issues if these are NULL.
    dispersion = suppressWarnings(mod.fit$rho$dispersion),
    theta = suppressWarnings(mod.fit$rho$theta),
    precision = suppressWarnings(mod.fit$rho$precision)
  )

  # 2. Use the returned linkinv function to calculate mu
  mu <- f_link$linkinv(eta)

  # =========================================================================

  tgt[[fitted_name]] <- mu

  # write back
  if (store_in == "model") {
    mod.fit$model <- tgt
  } else {
    mod.fit <- tgt
  }

  mod.fit
}
