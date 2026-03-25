#' Internal: get variance / mu' / linkinv for a family
#' @keywords internal
gee_family_fns <- function(family,
                           link = NULL,
                           dispersion = NULL,
                           theta = NULL,
                           precision = NULL,
                           zi_prob = NULL,
                           clamp_eps = 1e-8,
                           # optional overrides (advanced users)
                           varFn = NULL,
                           muprimeFn = NULL,     # mu.eta(eta)
                           linkinvFn = NULL,     # linkinv(eta)
                           ...) {

  # ----------------------------
  # Resolve family name + object
  # ----------------------------
  fam_obj <- NULL
  fam_name <- NULL

  if (inherits(family, "family")) {
    fam_obj <- family
    fam_name <- fam_obj$family
    if (is.null(link)) link <- fam_obj$link

  } else if (is.character(family) && length(family) == 1L) {

    fam_name <- family

    fam_fun <- tryCatch(get(family, mode = "function"), error = function(e) NULL)
    if (!is.null(fam_fun)) {
      fam_obj <- tryCatch({
        fmls <- names(formals(fam_fun))
        if (!is.null(link) && "link" %in% fmls) fam_fun(link = link) else fam_fun()
      }, error = function(e) NULL)

      if (!is.null(fam_obj) && inherits(fam_obj, "family")) {
        fam_name <- fam_obj$family
        if (is.null(link)) link <- fam_obj$link
      } else {
        fam_obj <- NULL
      }
    }

    if (is.null(link)) {
      fam_key <- tolower(gsub("[[:space:]_\\-\\.]", "", fam_name))
      link <- switch(fam_key,
                     gaussian = "identity",
                     normal   = "identity",
                     binomial = "logit",
                     quasibinomial = "logit",
                     poisson  = "log",
                     quasipoisson = "log",
                     gamma = "inverse",
                     inversegaussian = "1/mu^2",
                     "identity"
      )
    }

  } else {
    stop("family must be a 'family' object or a character scalar.")
  }

  if (!is.character(link) || length(link) != 1L) stop("link must be a scalar character.")
  link_lc <- tolower(link)

  clamp_eps <- as.numeric(clamp_eps)
  if (!is.finite(clamp_eps) || clamp_eps < 0) clamp_eps <- 0

  # ----------------------------
  # Link pieces (eta -> mu) and mu.eta
  # ----------------------------
  if (!is.null(linkinvFn)) {
    linkinv <- linkinvFn
  } else if (!is.null(fam_obj) && is.function(fam_obj$linkinv)) {
    linkinv <- fam_obj$linkinv
  } else {
    lk <- stats::make.link(link_lc)
    linkinv <- lk$linkinv
  }

  if (!is.null(muprimeFn)) {
    mu_eta <- muprimeFn
  } else if (!is.null(fam_obj) && is.function(fam_obj$mu.eta)) {
    mu_eta <- fam_obj$mu.eta
  } else {
    lk <- stats::make.link(link_lc)
    mu_eta <- lk$mu.eta
  }

  # ---- numeric safety wrappers for some links ----
  if (identical(link_lc, "log")) {
    linkinv <- function(eta) exp(pmin(eta, 700))
    mu_eta  <- function(eta) exp(pmin(eta, 700))
  }
  # if (identical(link_lc, "log")) {
  #
  #   logxmax <- log(.Machine$double.xmax)
  #
  #   # default: exp(eta) itself is safe
  #   eta_cap <- logxmax - 1
  #
  #   # families that later square/cube mu in variance:
  #   #   Gamma: Var ~ mu^2
  #   #   NB2:   Var ~ mu + mu^2/theta
  #   #   InvGauss: Var ~ mu^3
  #   if (fam_key %in% c("gamma","quasigamma") ||
  #       fam_key %in% c("nb","negbinomial","negativebinomial") ||
  #       grepl("negativebinomial", fam_key, fixed = TRUE)) {
  #     eta_cap <- (logxmax / 2) - 1
  #   }
  #   if (fam_key %in% c("inversegaussian","inversegaussianfamily")) {
  #     eta_cap <- (logxmax / 3) - 1
  #   }
  #
  #   linkinv <- function(eta) exp(pmin(eta, eta_cap))
  #   mu_eta  <- function(eta) exp(pmin(eta, eta_cap))
  # }

  if (identical(link_lc, "cloglog")) {
    linkinv <- function(eta) {
      e <- exp(pmin(eta, 700))
      1 - exp(-e)
    }
    mu_eta <- function(eta) {
      e <- exp(pmin(eta, 700))
      e * exp(-e)
    }
  }

  if (identical(link_lc, "inverse")) {
    tol_eta <- sqrt(.Machine$double.eps)
    linkinv <- function(eta) {
      eta <- as.numeric(eta)
      eta_safe <- ifelse(abs(eta) < tol_eta, ifelse(eta >= 0, tol_eta, -tol_eta), eta)
      1 / eta_safe
    }
    mu_eta <- function(eta) {
      eta <- as.numeric(eta)
      eta_safe <- ifelse(abs(eta) < tol_eta, ifelse(eta >= 0, tol_eta, -tol_eta), eta)
      -1 / (eta_safe^2)
    }
  }

  if (identical(link_lc, "1/mu^2")) {
    tol_eta <- sqrt(.Machine$double.eps)
    linkinv <- function(eta) {
      eta <- pmax(as.numeric(eta), tol_eta)
      1 / sqrt(eta)
    }
    mu_eta <- function(eta) {
      eta <- pmax(as.numeric(eta), tol_eta)
      -(0.5) * eta^(-3/2)
    }
  }

  # ---- IMPORTANT: preserve matrix dimensions for CV ----
  # CV passes eta as an nobs x G matrix. Many base ops (ifelse/as.numeric)
  # drop dim(). This wrapper guarantees linkinv/mu_eta return same shape.
  .wrap_dim <- function(fun) {
    force(fun)
    function(x) {
      d <- dim(x)
      out <- fun(as.vector(x))     # always apply elementwise on vector
      out <- as.numeric(out)
      if (!is.null(d)) dim(out) <- d
      out
    }
  }

  linkinv <- .wrap_dim(linkinv)
  mu_eta  <- .wrap_dim(mu_eta)


  # ----------------------------
  # Nuisance-aware variance varFn(mu)
  # ----------------------------
  fam_key <- tolower(gsub("[[:space:]_\\-\\.]", "", fam_name))

  is_nb <- (fam_key %in% c("nb","negbinomial","negativebinomial") ||
              grepl("negativebinomial", fam_key, fixed = TRUE))

  is_beta <- (fam_key %in% c("beta", "betar","betaregression") ||
                grepl("betareg", fam_key, fixed = TRUE))

  if (!is.null(varFn)) {

    vfun <- varFn

  } else {

    if (!is.null(fam_obj) && is.function(fam_obj$variance) && !is_nb && !is_beta) {
      var0 <- fam_obj$variance
    } else {
      var0 <- switch(fam_key,
                     gaussian = function(mu) rep.int(1, length(mu)),
                     normal   = function(mu) rep.int(1, length(mu)),
                     binomial = function(mu) mu * (1 - mu),
                     quasibinomial = function(mu) mu * (1 - mu),
                     poisson  = function(mu) mu,
                     quasipoisson = function(mu) mu,
                     gamma = function(mu) mu^2,
                     inversegaussian = function(mu) mu^3,
                     inversegaussianfamily = function(mu) mu^3,
                     function(mu) stop("No built-in variance for family='", fam_name,
                                       "'. Supply varFn=, or pass a family object with a $variance function.")
      )
    }

    vfun <- function(mu) {

      mu <- as.numeric(mu)
      if (any(!is.finite(mu))) stop("Non-finite mu encountered in varFn().")

      if (fam_key %in% c("binomial","quasibinomial") || is_beta) {
        if (clamp_eps > 0) mu <- pmin(pmax(mu, clamp_eps), 1 - clamp_eps)
      } else {
        if (clamp_eps > 0) mu <- pmax(mu, clamp_eps)
      }

      if (is_beta) {
        if (is.null(precision)) stop("beta regression variance requires 'precision' (phi).")
        ph <- as.numeric(precision)
        if (!is.finite(ph) || ph <= 0) stop("precision must be positive/finite.")
        return((mu * (1 - mu)) / (1 + ph))
      }

      if (is_nb) {
        if (!is.null(theta)) {
          th <- as.numeric(theta)
          if (!is.finite(th) || th <= 0) stop("theta must be positive/finite for NB variance.")
          return(mu + (mu^2) / th)
        }
        if (!is.null(fam_obj) && is.function(fam_obj$variance)) {
          return(as.numeric(fam_obj$variance(mu)))
        }
        stop("Negative binomial variance requires theta= or a NB family object.")
      }

      vv <- as.numeric(var0(mu))

      if (!is.null(dispersion)) {
        ph <- as.numeric(dispersion)
        if (!is.finite(ph) || ph <= 0) stop("dispersion must be positive/finite.")
        vv <- ph * vv
      }

      if (any(!is.finite(vv)) || any(vv <= 0)) stop("Non-finite or non-positive variance returned.")
      vv
    }
  }

  out <- list(
    family  = fam_name,
    link    = link_lc,
    linkinv = linkinv,
    mu_eta  = mu_eta,
    varFn   = vfun,
    variance = vfun
  )

  out$linkinvFn  <- out$linkinv
  out$muprimeFn  <- out$mu_eta

  out
}


# ============================================================================
# HELPER for Family and Nuisance Parameters
# ============================================================================
#' @keywords internal
#' @noRd
get_family_info <- function(fit.initial) {

  family_obj   <- fit.initial$family
  fam_name_raw <- family_obj$family

  # ---- helpers ----
  .num1 <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) < 1L) return(NA_real_)
    x[1]
  }
  .is_pos <- function(x) is.finite(x) && (x > 0)
  .key <- function(x) tolower(gsub("[^a-z0-9]", "", x))

  fam_key_raw   <- tolower(fam_name_raw)
  fam_key_clean <- .key(fam_name_raw)

  # link (needed for C++; harmless for R)
  link_cpp <- tolower(if (!is.null(family_obj$link)) family_obj$link else fit.initial$family$link)

  # --- family key -> C++ family string (NULL => C++ disabled, R still works)
  family_cpp <- NULL

  if (startsWith(fam_key_raw, "negative binomial") ||
      grepl("negativebinomial", fam_key_clean, fixed = TRUE) ||
      grepl("negbinomial", fam_key_clean, fixed = TRUE) ||
      grepl("negbin", fam_key_clean, fixed = TRUE) ||
      identical(fam_key_clean, "nb")) {
    family_cpp <- "negbinomial"

  } else if (fam_key_clean %in% c("gaussian", "normal")) {
    family_cpp <- "gaussian"

  } else if (fam_key_clean %in% c("binomial", "quasibinomial")) {
    family_cpp <- fam_key_clean

  } else if (fam_key_clean %in% c("poisson", "quasipoisson")) {
    family_cpp <- fam_key_clean

  } else if (startsWith(fam_key_raw, "Gamma") ||
             startsWith(fam_key_raw, "gamma") ||
             startsWith(fam_key_raw, "quasigamma") ){
  #(fam_key_clean %in% c("gamma", "quasigamma")) {
    family_cpp <- "gamma"
  } else if (startsWith(fam_key_raw, "beta")) {
# fam_key_clean %in% c("beta", "betar", "betaregression") ||
             # grepl("betareg", fam_key_clean, fixed = TRUE))
    family_cpp <- "beta"
  }

  # default nuisance scalar
  dispersion_cpp <- 1.0

  # -------------------------
  # nuisance extraction
  # -------------------------
  if (identical(family_cpp, "negbinomial")) {

    th0 <- NA_real_

    # 1) preferred: getTheta(trans=TRUE) or getTheta(TRUE)
    if (is.function(family_obj$getTheta)) {
      th0 <- tryCatch({
        fmls <- names(formals(family_obj$getTheta))
        if ("trans" %in% fmls) family_obj$getTheta(trans = TRUE)
        else family_obj$getTheta(TRUE)
      }, error = function(e) NA_real_)
      th0 <- .num1(th0)
    }

    # 2) fallback slots
    if (!.is_pos(th0)) th0 <- .num1(family_obj$theta)
    if (!.is_pos(th0)) th0 <- .num1(fit.initial$family$theta)

    # 3) parse from family string "Negative Binomial(3.2)"
    if (!.is_pos(th0)) {
      th1 <- .num1(sub(".*\\(([^\\)]+)\\).*", "\\1", fam_name_raw))
      if (.is_pos(th1)) th0 <- th1
    }

    if (.is_pos(th0)) dispersion_cpp <- th0 else dispersion_cpp <- 1.0

  } else if (identical(family_cpp, "beta")) {

    ph0 <- NA_real_

    # common fit slots
    ph0 <- .num1(fit.initial$precision)
    if (!.is_pos(ph0)) ph0 <- .num1(fit.initial$phi)

    # family slots (custom families sometimes store it here)
    if (!.is_pos(ph0)) ph0 <- .num1(family_obj$precision)
    if (!.is_pos(ph0)) ph0 <- .num1(family_obj$phi)
    if (!.is_pos(ph0)) ph0 <- .num1(family_obj$theta)

    # some families provide getTheta()
    if (!.is_pos(ph0) && is.function(family_obj$getTheta)) {
      ph0 <- tryCatch({
        fmls <- names(formals(family_obj$getTheta))
        if ("trans" %in% fmls) family_obj$getTheta(trans = TRUE)
        else family_obj$getTheta(TRUE)
      }, error = function(e) NA_real_)
      ph0 <- .num1(ph0)
    }

    if (.is_pos(ph0)) dispersion_cpp <- ph0 else dispersion_cpp <- 1.0
  }

  list(
    family_cpp      = family_cpp,      # NULL => C++ disabled, R still runs
    dispersion_cpp  = dispersion_cpp,  # theta for NB, phi for beta, else 1.0
    link_cpp        = link_cpp
  )
}


#' Internal: compute eta = X %*% beta with optional "accumulate" to avoid dense X
#' @keywords internal
.dt_linpred <- function(dt, xcols, beta,
                        out = "eta",
                        method = c("accumulate", "matrix"),
                        return = c("dt", "vector")) {

  method <- match.arg(method)
  return <- match.arg(return)

  dt <- data.table::as.data.table(dt)

  xcols <- as.character(xcols)
  beta  <- as.numeric(beta)

  if (length(xcols) != length(beta)) {
    stop("length(xcols) must equal length(beta).")
  }
  miss <- setdiff(xcols, names(dt))
  if (length(miss)) stop("Missing xcols in dt: ", paste(miss, collapse = ", "))

  n <- nrow(dt)
  if (n == 0L) {
    if (return == "vector") return(numeric(0))
    dt[, (out) := numeric(0)]
    return(dt)
  }

  if (method == "matrix") {
    X <- as.matrix(dt[, ..xcols])
    eta <- as.numeric(X %*% beta)
  } else {
    eta <- numeric(n)
    for (j in seq_along(xcols)) {
      bj <- beta[j]
      if (!is.finite(bj) || bj == 0) next
      eta <- eta + dt[[xcols[j]]] * bj
    }
  }

  if (return == "vector") return(eta)

  dt[, (out) := eta]
  dt
}


# -----------------------------------------------------------------------------
# Core updater: compute eta, p (=mu), muprime, v, sqrtv, resid (scaled by 1/sd)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Core updater: eta, mu (=p), muprime, v, sqrtv, resid
#  - uses .dt_linpred()
#  - supports family as a 'family' object OR a family name string
#  - supports linpred_method = c("accumulate","matrix")
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
fgee_update_working_cols_dt <- function(dx,
                                        namesd,
                                        beta,
                                        family,
                                        link = NULL,
                                        exact = FALSE,
                                        gaussian_v_by = NULL,
                                        update_nuisance = c("fixed","moment"),
                                        dispersion = NULL,
                                        theta = NULL,
                                        precision = NULL,
                                        zi_prob = NULL,
                                        clamp_eps = 1e-8,
                                        # linear predictor control
                                        linpred_method = c("accumulate", "matrix"),
                                        # optional blockwise eta
                                        eta_by = NULL,
                                        # output column names
                                        eta_col = "eta",
                                        p_col = "p",
                                        v_col = "v",
                                        sqrtv_col = "sqrtv",
                                        muprime_col = "muprime",
                                        resid_col = "resid",
                                        copy = TRUE) {

  update_nuisance <- match.arg(update_nuisance)
  linpred_method  <- match.arg(linpred_method)

  dd <- data.table::as.data.table(dx)
  if (copy) dd <- data.table::copy(dd)

  namesd <- as.character(namesd)
  beta   <- as.numeric(beta)

  if (length(beta) != length(namesd)) {
    stop("beta length (", length(beta), ") must equal length(namesd) (", length(namesd), ").")
  }
  if (!("Y" %in% names(dd))) stop("dx must contain column 'Y'.")

  # --- family name string (for nuisance updater + clamping decisions)
  fam_name <- if (inherits(family, "family")) family$family else as.character(family)
 # fam_key  <- tolower(gsub("[[:space:]_\\-\\.]", "", fam_name))
  fam_key <- tolower(gsub("[^a-z]", "", tolower(fam_name))) # would need to update this if we used tweeide as it may not properly strip everything off

  is_nb <- (fam_key %in% c("nb","negbinomial","negativebinomial") ||
              grepl("negativebinomial", fam_key, fixed = TRUE))
  is_beta <- (fam_key %in% c("beta", "betar","betaregression") ||
                grepl("betareg", fam_key, fixed = TRUE))

  # --- build link pieces (eta -> mu) + mu.eta; nuisance values don't matter here
  f_link <- gee_family_fns(
    family = family,
    link = link,
    dispersion = dispersion,
    theta = theta,
    precision = precision,
    zi_prob = zi_prob,
    clamp_eps = clamp_eps
  )

  # -----------------------
  # eta (optionally blockwise)
  # -----------------------
  if (is.null(eta_by)) {
    dd[, (eta_col) := .dt_linpred(dd, xcols = namesd, beta = beta,
                                  method = linpred_method, return = "vector")]
  } else {
    dd[, (eta_col) := .dt_linpred(.SD, xcols = namesd, beta = beta,
                                  method = linpred_method, return = "vector"),
       by = eta_by, .SDcols = namesd]
  }

  eta <- dd[[eta_col]]

  # -----------------------
  # mu = linkinv(eta)
  # -----------------------
  mu <- as.numeric(f_link$linkinv(eta))

  # Clamp mu where needed (keeps varFn + muprime stable)
  if (clamp_eps > 0) {
    if (fam_key %in% c("binomial","quasibinomial") || is_beta) {
      mu <- pmin(pmax(mu, clamp_eps), 1 - clamp_eps)
    } else if (fam_key %in% c("poisson","quasipoisson","gamma","inversegaussian","inversegaussianfamily") || is_nb) {
      mu <- pmax(mu, clamp_eps)
    }
  }

  dd[, (p_col) := mu]

  # -----------------------
  # optional nuisance update using UPDATED mu
  # -----------------------
  if (update_nuisance == "moment") {
    nu <- .fgee_update_nuisance_moments(
      y = dd$Y, mu = mu,
      family = fam_name,
      dispersion = dispersion,
      theta = theta,
      precision = precision,
      clamp_eps = clamp_eps
    )
    dispersion <- nu$dispersion
    theta      <- nu$theta
    precision  <- nu$precision
  }

  # -----------------------
  # build nuisance-aware variance function varFn(mu)
  # -----------------------
  f <- gee_family_fns(
    family = family,
    link = f_link$link,   # ensure consistent with resolved link
    dispersion = dispersion,
    theta = theta,
    precision = precision,
    zi_prob = zi_prob,
    clamp_eps = clamp_eps
  )

  # -----------------------
  # v
  # -----------------------
  if (tolower(f$family) %in% c("gaussian","normal") && !is.null(gaussian_v_by)) {

    if (!(gaussian_v_by %in% names(dd))) stop("gaussian_v_by not found: ", gaussian_v_by)

    dd[, .res_raw_tmp := (Y - get(p_col))]
    dd[, (v_col) := stats::var(.res_raw_tmp), by = gaussian_v_by]
    dd[, .res_raw_tmp := NULL]

    # fallback for singleton / NA groups
    if (any(!is.finite(dd[[v_col]]))) {
      v_global <- stats::var(dd$Y - dd[[p_col]], na.rm = TRUE)
      dd[!is.finite(get(v_col)), (v_col) := v_global]
    }
    if (any(dd[[v_col]] <= 0, na.rm = TRUE)) stop("Non-positive gaussian variance encountered.")

  } else {

    vv <- as.numeric(f$varFn(mu))
    if (any(!is.finite(vv)) || any(vv <= 0)) stop("Bad v from varFn; check parameters and mu range.")
    dd[, (v_col) := vv]
  }

  # -----------------------
  # sqrtv (WRITE IT deterministically)
  # -----------------------
  dd[, (sqrtv_col) := sqrt(get(v_col))]

  # -----------------------
  # muprime = d mu / d eta
  # -----------------------
  muprime <- as.numeric(f_link$mu_eta(eta))

  # If link is one we can express stably in terms of mu, do it (keeps clamp consistent)
  link_lc <- tolower(f_link$link)
  if (link_lc %in% c("identity","log","logit","probit","cloglog","inverse","sqrt","1/mu^2")) {
    muprime2 <- try(.fgee_muprime_from_mu(mu, link = link_lc, clamp_eps = clamp_eps), silent = TRUE)
    if (!inherits(muprime2, "try-error")) muprime <- muprime2
  }

  dd[, (muprime_col) := as.numeric(muprime)]


  # -----------------------
  # resid standardized by sqrt(v)
  # -----------------------
  if (isTRUE(exact)) {
    dd[, (resid_col) := Y / get(sqrtv_col)]
  } else {
    dd[, (resid_col) := (Y - get(p_col)) / get(sqrtv_col)]
  }

  # optional: attach nuisance info (handy for debugging; harmless)
  attr(dd, "nuisance") <- list(dispersion = dispersion, theta = theta, precision = precision, zi_prob = zi_prob)

  dd
}


#' @keywords internal
#' @noRd
.fgee_muprime_from_mu <- function(mu, link, clamp_eps = 1e-8) {
  link <- tolower(link)
  mu <- as.numeric(mu)

  if (!is.finite(clamp_eps) || clamp_eps < 0) clamp_eps <- 0

  clamp01 <- function(x) {
    if (clamp_eps <= 0) return(x)
    pmin(pmax(x, clamp_eps), 1 - clamp_eps)
  }

  switch(link,
         "identity" = rep.int(1, length(mu)),

         "log" = pmax(mu, 0),

         "logit" = {
           m <- clamp01(mu)
           m * (1 - m)
         },

         "probit" = {
           m <- clamp01(mu)
           eta <- stats::qnorm(m)
           stats::dnorm(eta)
         },

         "cloglog" = {
           m <- clamp01(mu)
           om <- pmax(1 - m, max(clamp_eps, 1e-15))  # protect log(0)
           om * (-log(om))                          # (1-mu) * exp(eta) where exp(eta)=-log(1-mu)
         },

         "inverse" = {
           # mu = 1/eta  => dmu/deta = -1/eta^2 = -mu^2
           -(mu^2)
         },

         # not requested, but useful + already supported by make.link:
         "sqrt" = 2 * sqrt(pmax(mu, 0)),

         "1/mu^2" = {
           # eta = 1/mu^2 => mu = eta^{-1/2} => dmu/deta = -(1/2)*eta^{-3/2} = -(1/2)*mu^3 (mu>0)
           -(0.5) * (pmax(mu, 0)^3)
         },

         stop("Unsupported link for .fgee_muprime_from_mu(): ", link)
  )
}

#
# # ---- helpers ----
# .fgee_safe_exp <- function(x) {
#   # avoid overflow; exp(709) ~ 8e307 near double max
#   exp(pmin(x, 700))
# }
#
#
# .fgee_linkinv <- function(eta, link) {
#   link <- tolower(link)
#   switch(link,
#          "identity" = eta,
#          "log"      = exp(pmin(eta, 700)),  # avoid overflow
#          "logit"    = stats::plogis(eta),
#          stop("Unsupported link: ", link, ". Only identity/log/logit are supported.")
#   )
# }

# ------------------------------------------------------------
# Optional: moment updates for nuisance params using *current* mu
# (this is the piece you were worried about after 1-step beta changes)
# ------------------------------------------------------------
#' @keywords internal
#' @noRd
.fgee_update_nuisance_moments <- function(y, mu, family,
                                          dispersion = NULL,
                                          theta = NULL,
                                          precision = NULL,
                                          clamp_eps = 1e-8) {

  fam0 <- tolower(gsub("[^a-z]", "", tolower(family)))

  is_nb   <- identical(fam0, "nb") ||
    grepl("negbinomial", fam0, fixed = TRUE) ||
    grepl("negativebinomial", fam0, fixed = TRUE)

  is_beta <- (fam0 %in% c("beta", "betar", "betaregression")) ||
    grepl("betareg", fam0, fixed = TRUE)

  is_quasi <- fam0 %in% c("quasipoisson", "quasibinomial")

  out <- list(dispersion = dispersion, theta = theta, precision = precision)

  # quasi dispersion (Pearson-ish)
  if (is_quasi) {
    V0 <- if (fam0 == "quasipoisson") pmax(mu, clamp_eps) else {
      m <- pmin(pmax(mu, clamp_eps), 1 - clamp_eps)
      m * (1 - m)
    }
    phi_hat <- mean(((y - mu)^2) / pmax(V0, clamp_eps), na.rm = TRUE)
    if (is.finite(phi_hat) && phi_hat > 0) out$dispersion <- phi_hat
  }

  # Negative binomial: Var = mu + mu^2/theta = mu + alpha*mu^2, alpha=1/theta
  if (is_nb) {
    denom <- sum(mu^2, na.rm = TRUE)
    if (is.finite(denom) && denom > 0) {
      alpha_hat <- sum(((y - mu)^2 - y), na.rm = TRUE) / denom
      alpha_hat <- max(alpha_hat, 0)
      if (alpha_hat > 0) out$theta <- 1 / alpha_hat
    }
  }

  # Beta regression: Var = mu(1-mu)/(1+phi)
  # if (is_beta) {
  #   m <- pmin(pmax(mu, clamp_eps), 1 - clamp_eps)
  #   denom <- pmax(m * (1 - m), clamp_eps)
  #   inv1p <- mean(((y - m)^2) / denom, na.rm = TRUE)
  #   if (is.finite(inv1p) && inv1p > 0) {
  #     out$precision <- max(clamp_eps, (1 / inv1p) - 1)
  #   }
  # }

  # Beta regression: Var = mu(1-mu)/(1+phi)
  if (is_beta) {
    m <- pmin(pmax(mu, clamp_eps), 1 - clamp_eps)
    denom <- pmax(m * (1 - m), clamp_eps)
    inv1p <- mean(((y - m)^2) / denom, na.rm = TRUE)

    if (is.finite(inv1p) && inv1p > 0) {
      # floor away from 0 so phi can't be Inf
      inv1p <- max(inv1p, clamp_eps)
      phi_hat <- (1 / inv1p) - 1
      if (is.finite(phi_hat) && phi_hat > 0) {
        out$precision <- max(clamp_eps, phi_hat)
      }
    }
  }


  out
}


