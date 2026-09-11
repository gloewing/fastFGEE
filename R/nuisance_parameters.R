# Safe scalar nuisance-parameter handling ---------------------------------
#
# The mean-model coefficient update remains a one-step GEE update. Scalar
# variance-function nuisance parameters are extracted/estimated once for the
# initial working covariance, held fixed throughout smoothing-parameter
# selection, and updated once at the final coefficient estimate for inference.
# They are never re-estimated for every lambda candidate.

.fgee_nuisance_or <- function(x, y) if (is.null(x)) y else x

.fgee_nuisance_family_key <- function(family) {
  fam <- if (is.list(family) && !is.null(family$family)) {
    family$family
  } else {
    as.character(family)[1L]
  }
  z <- tolower(trimws(fam))
  if (grepl("negative[ -]?binomial|negbin|nbinom", z)) return("negbinomial")
  if (grepl("beta", z) && !grepl("binomial", z)) return("beta")
  if (grepl("gamma", z)) return("gamma")
  "none"
}

.fgee_positive_scalar <- function(x) {
  is.numeric(x) && length(x) >= 1L && is.finite(x[1L]) && x[1L] > 0
}

# Precondition check for the families whose working variance is undefined
# without a scalar nuisance. Returns the resolved value invisibly so callers can
# record it; errors rather than letting a silent 1.0 default through.
#
#' @keywords internal
#' @noRd
.fgee_assert_nuisance_frozen <- function(glmfit, nuisance = NULL) {
  key <- .fgee_nuisance_family_key(glmfit$family)
  param <- .fgee_nuisance_parameter_name(key)
  if (is.null(param) || !key %in% c("negbinomial", "beta")) {
    return(invisible(NULL))
  }
  val <- .fgee_nuisance_value(nuisance, param, NULL)
  if (!.fgee_positive_scalar(val)) {
    val <- .fgee_nuisance_value(glmfit, param, NULL)
  }
  if (!.fgee_positive_scalar(val)) {
    stop(
      "A positive ", param, " is required to form the working variance for ",
      glmfit$family$family, "/", glmfit$family$link,
      ", and none could be resolved. Fit with mgcv::nb(theta = ) or ",
      "mgcv::betar(theta = ) so the value is carried on the family object.",
      call. = FALSE
    )
  }
  invisible(as.numeric(val)[1L])
}

.fgee_nuisance_parameter_name <- function(key) {
  switch(key,
    gamma = "dispersion",
    beta = "precision",
    negbinomial = "theta",
    NULL
  )
}

.fgee_nuisance_bounds <- function(key, control = list()) {
  defaults <- switch(key,
    gamma = c(1e-8, 1e4),
    beta = c(1e-3, 1e8),
    negbinomial = c(1e-4, 1e8),
    stop("Unsupported nuisance family.")
  )
  b <- .fgee_nuisance_or(control$bounds, defaults)
  b <- as.numeric(b)
  if (length(b) != 2L || any(!is.finite(b)) || b[1L] <= 0 || b[2L] <= b[1L]) {
    stop("Nuisance bounds must be two finite positive increasing values.")
  }
  b
}

.fgee_nuisance_weights <- function(n, cluster = NULL, weights = NULL,
                                   weighting = c("observation", "cluster_equal")) {
  weighting <- match.arg(weighting)
  if (!is.null(weights)) {
    w <- as.numeric(weights)
    if (length(w) != n || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
      stop("Nuisance weights must be finite, non-negative, and match y.")
    }
  } else if (weighting == "cluster_equal" && !is.null(cluster)) {
    cluster <- as.character(cluster)
    if (length(cluster) != n || anyNA(cluster)) {
      stop("cluster must match y and contain no missing values.")
    }
    tab <- table(cluster)
    w <- 1 / as.numeric(tab[cluster])
  } else {
    w <- rep.int(1, n)
  }
  # Normalize for numerically comparable objective values. The minimizer is
  # unchanged by this constant rescaling.
  w * (n / sum(w))
}

.fgee_prepare_nuisance_data <- function(y, mu, key, cluster = NULL,
                                        weights = NULL, control = list()) {
  y <- as.numeric(y)
  mu <- as.numeric(mu)
  if (length(y) != length(mu) || !length(y)) {
    stop("y and mu must be non-empty and have the same length.")
  }
  keep <- is.finite(y) & is.finite(mu)
  if (!all(keep)) {
    warning(sum(!keep), " non-finite y/mu rows were omitted from nuisance estimation.",
            call. = FALSE)
    y <- y[keep]
    mu <- mu[keep]
    if (!is.null(cluster)) cluster <- cluster[keep]
    if (!is.null(weights)) weights <- weights[keep]
  }
  eps <- as.numeric(.fgee_nuisance_or(control$eps, 1e-8))[1L]
  if (!is.finite(eps) || eps <= 0 || eps >= 0.1) stop("control$eps is invalid.")
  clipped_y <- 0L
  clipped_mu <- 0L

  if (key == "gamma") {
    if (any(y <= 0)) stop("Gamma nuisance estimation requires strictly positive outcomes.")
    clipped_mu <- sum(mu <= eps)
    mu <- pmax(mu, eps)
  } else if (key == "beta") {
    action <- .fgee_nuisance_or(control$beta_boundary, "warn_clip")
    bad <- y <= 0 | y >= 1
    clipped_y <- sum(bad)
    if (clipped_y) {
      if (identical(action, "error")) {
        stop("Beta nuisance estimation requires outcomes strictly inside (0, 1).")
      }
      if (!identical(action, "clip") && !identical(action, "warn_clip")) {
        stop("control$beta_boundary must be 'error', 'clip', or 'warn_clip'.")
      }
      if (identical(action, "warn_clip")) {
        warning(clipped_y, " beta outcomes were clipped into (eps, 1-eps); ",
                "this adjustment is recorded in the nuisance diagnostics.", call. = FALSE)
      }
      y <- pmin(pmax(y, eps), 1 - eps)
    }
    clipped_mu <- sum(mu <= eps | mu >= 1 - eps)
    mu <- pmin(pmax(mu, eps), 1 - eps)
  } else if (key == "negbinomial") {
    if (any(y < 0) || any(abs(y - round(y)) > 1e-7)) {
      stop("Negative-binomial nuisance estimation requires non-negative integer outcomes.")
    }
    clipped_mu <- sum(mu <= eps)
    mu <- pmax(mu, eps)
  }

  weighting <- .fgee_nuisance_or(control$weighting, "observation")
  w <- .fgee_nuisance_weights(length(y), cluster, weights, weighting)
  list(y = y, mu = mu, w = w, cluster = cluster,
       clipped_y = clipped_y, clipped_mu = clipped_mu, eps = eps,
       weighting = weighting)
}

.fgee_nuisance_moment <- function(dat, key, bounds) {
  y <- dat$y; mu <- dat$mu; w <- dat$w
  sw <- sum(w)
  if (key == "gamma") {
    value <- sum(w * ((y - mu) / mu)^2) / sw
    boundary <- FALSE
  } else if (key == "beta") {
    # Solves sum{(y-mu)^2 - mu(1-mu)/(1+precision)} = 0. This is
    # more stable near the boundaries than averaging row-wise ratios.
    num <- sum(w * mu * (1 - mu))
    den <- sum(w * (y - mu)^2)
    value <- num / max(den, .Machine$double.eps) - 1
    boundary <- !is.finite(value) || value <= bounds[1L] || value >= bounds[2L]
  } else if (key == "negbinomial") {
    # NB2: Var(Y)=mu+mu^2/theta. A non-positive over-dispersion
    # estimate is the Poisson boundary rather than a reason to silently reuse
    # an unrelated previous theta.
    num <- sum(w * mu^2)
    den <- sum(w * ((y - mu)^2 - mu))
    if (!is.finite(den) || den <= 0) {
      value <- bounds[2L]
      boundary <- TRUE
    } else {
      value <- num / den
      boundary <- value <= bounds[1L] || value >= bounds[2L]
    }
  } else stop("Unsupported nuisance family.")

  if (!is.finite(value) || value <= 0) {
    return(list(ok = FALSE, value = NA_real_, boundary = TRUE,
                message = "Moment nuisance estimate was not finite and positive."))
  }
  value <- min(max(value, bounds[1L]), bounds[2L])
  list(ok = TRUE, value = value, boundary = boundary, message = NULL)
}

.fgee_nuisance_profile_objective <- function(log_value, dat, key) {
  value <- exp(log_value)
  y <- dat$y; mu <- dat$mu; w <- dat$w
  ll <- switch(key,
    gamma = stats::dgamma(y, shape = 1 / value, scale = value * mu, log = TRUE),
    beta = {
      a <- mu * value
      b <- (1 - mu) * value
      -lbeta(a, b) + (a - 1) * log(y) + (b - 1) * log1p(-y)
    },
    negbinomial = stats::dnbinom(y, size = value, mu = mu, log = TRUE),
    stop("Unsupported nuisance family.")
  )
  if (any(!is.finite(ll))) return(.Machine$double.xmax / 100)
  -sum(w * ll) / sum(w)
}

.fgee_nuisance_profile <- function(dat, key, bounds, control = list()) {
  tol <- as.numeric(.fgee_nuisance_or(control$tol, 1e-7))[1L]
  opt <- tryCatch(
    stats::optimize(
      f = .fgee_nuisance_profile_objective,
      interval = log(bounds),
      dat = dat,
      key = key,
      tol = tol
    ),
    error = identity
  )
  if (inherits(opt, "error") || !is.finite(opt$minimum) || !is.finite(opt$objective)) {
    return(list(ok = FALSE, value = NA_real_, boundary = FALSE,
                objective = NA_real_, message = if (inherits(opt, "error")) conditionMessage(opt) else "Profile optimization failed."))
  }
  value <- exp(opt$minimum)
  span <- diff(log(bounds))
  boundary <- (opt$minimum - log(bounds[1L]) <= max(1e-6, span * 1e-4)) ||
              (log(bounds[2L]) - opt$minimum <= max(1e-6, span * 1e-4))
  list(ok = TRUE, value = value, boundary = boundary,
       objective = opt$objective, message = NULL)
}

# Estimate one scalar variance-function nuisance parameter while holding the
# fitted mean fixed. Profile likelihood here is marginal/composite under
# within-cluster dependence; the cluster weighting policy is explicit.
fgee_estimate_nuisance <- function(y, mu, family,
                                   method = c("profile", "moment"),
                                   cluster = NULL, weights = NULL,
                                   start = NULL, control = list()) {
  method <- match.arg(method)
  key <- .fgee_nuisance_family_key(family)
  if (key == "none") stop("This family has no supported scalar nuisance parameter.")
  dat <- .fgee_prepare_nuisance_data(y, mu, key, cluster, weights, control)
  bounds <- .fgee_nuisance_bounds(key, control)
  est <- if (method == "profile") {
    .fgee_nuisance_profile(dat, key, bounds, control)
  } else {
    .fgee_nuisance_moment(dat, key, bounds)
  }
  fallback <- FALSE
  if (!isTRUE(est$ok) && method == "profile") {
    est <- .fgee_nuisance_moment(dat, key, bounds)
    fallback <- TRUE
  }
  if (!isTRUE(est$ok) && .fgee_positive_scalar(start)) {
    est <- list(ok = TRUE, value = as.numeric(start)[1L], boundary = FALSE,
                objective = NA_real_, message = "Kept the previous valid nuisance value.")
    fallback <- TRUE
  }
  if (!isTRUE(est$ok)) {
    stop("Unable to obtain a finite positive nuisance estimate for ", key,
         ". Supply a valid initial value or inspect the mean fit.")
  }
  param <- .fgee_nuisance_parameter_name(key)
  value <- min(max(as.numeric(est$value)[1L], bounds[1L]), bounds[2L])
  ans <- list(
    family = key,
    parameter = param,
    value = value,
    method = if (fallback) paste0(method, "_fallback") else method,
    converged = TRUE,
    boundary = isTRUE(est$boundary),
    bounds = bounds,
    objective = .fgee_nuisance_or(est$objective, NA_real_),
    weighting = dat$weighting,
    n_used = length(dat$y),
    clipped_y = dat$clipped_y,
    clipped_mu = dat$clipped_mu,
    message = est$message
  )
  ans[[param]] <- value
  class(ans) <- c("fgee_nuisance", "list")
  ans
}

.fgee_nuisance_from_family <- function(family, key = .fgee_nuisance_family_key(family)) {
  val <- NULL
  source <- NULL
  if (is.list(family)) {
    if (key %in% c("beta", "negbinomial") && is.function(family$getTheta)) {
      val <- tryCatch(family$getTheta(TRUE), error = function(e) NULL)
      if (.fgee_positive_scalar(val)) source <- "family$getTheta(TRUE)"
    }
    param <- .fgee_nuisance_parameter_name(key)
    if (!.fgee_positive_scalar(val) && !is.null(param) && .fgee_positive_scalar(family[[param]])) {
      val <- family[[param]]
      source <- paste0("family$", param)
    }
  }
  if (!.fgee_positive_scalar(val)) return(NULL)
  param <- .fgee_nuisance_parameter_name(key)
  ans <- list(family = key, parameter = param, value = as.numeric(val)[1L],
              method = "model", source = source, converged = TRUE,
              boundary = FALSE, message = NULL)
  ans[[param]] <- ans$value
  class(ans) <- c("fgee_nuisance", "list")
  ans
}

.fgee_nuisance_from_fit <- function(fit) {
  if (is.null(fit)) return(NULL)
  fam <- .fgee_nuisance_or(fit$family,
                           if (!is.null(fit$model)) fit$model$family else NULL)
  key <- .fgee_nuisance_family_key(fam)
  if (key == "none") return(NULL)
  if (key == "gamma") {
    candidates <- list(fit$sig2, fit$scale,
                       if (is.list(fam)) fam$dispersion else NULL,
                       if (is.list(fam)) fam$scale else NULL)
    for (j in seq_along(candidates)) {
      if (.fgee_positive_scalar(candidates[[j]])) {
        val <- as.numeric(candidates[[j]])[1L]
        ans <- list(family = key, parameter = "dispersion", value = val,
                    dispersion = val, method = "model", source = c("fit$sig2", "fit$scale", "family$dispersion", "family$scale")[j],
                    converged = TRUE, boundary = FALSE, message = NULL)
        class(ans) <- c("fgee_nuisance", "list")
        return(ans)
      }
    }
    return(NULL)
  }
  .fgee_nuisance_from_family(fam, key)
}

.fgee_nuisance_legacy_values <- function(x) {
  if (is.null(x)) return(list())
  list(
    dispersion = x$dispersion,
    precision = x$precision,
    theta = x$theta,
    zi_prob = x$zi_prob
  )
}

.fgee_nuisance_value <- function(x, name, default = NULL) {
  if (is.null(x)) return(default)
  direct <- if (is.list(x)) x[[name]] else NULL
  if (.fgee_positive_scalar(direct)) return(as.numeric(direct)[1L])
  n <- attr(x, "nuisance", exact = TRUE)
  if (is.null(n) && is.list(x) && !is.null(x$nuisance)) n <- x$nuisance
  if (is.list(n)) {
    if (.fgee_positive_scalar(n[[name]])) return(as.numeric(n[[name]])[1L])
    if (identical(n$parameter, name) && .fgee_positive_scalar(n$value)) return(as.numeric(n$value)[1L])
    if (!is.null(n$final)) return(.fgee_nuisance_value(n$final, name, default))
  }
  default
}

.fgee_nuisance_variance <- function(mu, nuisance) {
  key <- nuisance$family
  value <- nuisance$value
  switch(key,
    gamma = value * mu^2,
    beta = mu * (1 - mu) / (1 + value),
    negbinomial = mu + mu^2 / value,
    stop("Unsupported nuisance family.")
  )
}

.fgee_apply_nuisance_to_working_data <- function(dd, nuisance,
                                                  y_col = "Y", mu_col = "p",
                                                  v_col = "v", sqrtv_col = "sqrtv",
                                                  resid_col = "resid") {
  if (!data.table::is.data.table(dd)) dd <- data.table::as.data.table(dd)
  if (!all(c(y_col, mu_col) %in% names(dd))) {
    stop("Working data must contain outcome and fitted-mean columns.")
  }
  mu <- as.numeric(dd[[mu_col]])
  v <- .fgee_nuisance_variance(mu, nuisance)
  if (any(!is.finite(v)) || any(v <= 0)) stop("Nuisance variance calculation was not positive and finite.")
  sqrtv <- sqrt(v)
  resid <- (as.numeric(dd[[y_col]]) - mu) / sqrtv
  data.table::set(dd, j = v_col, value = v)
  data.table::set(dd, j = sqrtv_col, value = sqrtv)
  data.table::set(dd, j = resid_col, value = resid)
  data.table::setattr(dd, "nuisance", nuisance)
  dd
}

.fgee_nuisance_map_dots <- function(core, dots) {
  nms <- names(dots)
  if (is.null(nms)) nms <- rep.int("", length(dots))
  fml <- names(formals(core))
  used <- nms[nzchar(nms)]
  next_pos <- 1L
  for (i in seq_along(dots)) {
    if (!nzchar(nms[i])) {
      while (next_pos <= length(fml) && fml[next_pos] %in% used) next_pos <- next_pos + 1L
      if (next_pos > length(fml) || identical(fml[next_pos], "...")) {
        stop("Unable to map an unnamed argument in nuisance integration wrapper.")
      }
      nms[i] <- fml[next_pos]
      used <- c(used, nms[i])
      next_pos <- next_pos + 1L
    }
  }
  names(dots) <- nms
  dots
}
