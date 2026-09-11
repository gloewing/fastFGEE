# Sandwich-scaled working restricted quasi-REML -------------------------------
#
# This is a smoothing-parameter criterion for the frozen one-step quadratic
# estimating equation.  It is not a likelihood for the original outcomes.

.fgee_qreml_start_grid <- function(q, coarse = c(-4, -2, 0, 2, 4),
                                    bound = 6, max_points = 625L) {
  q <- as.integer(q)
  coarse <- as.numeric(coarse)
  full_n <- length(coarse)^q
  if (full_n <= max_points) {
    out <- as.matrix(expand.grid(rep(list(coarse), q), KEEP.OUT.ATTRS = FALSE))
    return(out[apply(abs(out) <= bound, 1L, all), , drop = FALSE])
  }

  out <- matrix(0, 1L, q)
  for (a in coarse) out <- rbind(out, rep(a, q))
  for (j in seq_len(q)) {
    for (a in coarse[coarse != 0]) {
      u <- rep(0, q)
      u[j] <- a
      out <- rbind(out, u)
    }
  }
  unique(out[apply(abs(out) <= bound, 1L, all), , drop = FALSE])
}

#' Prepare working-H qREML from compact statistics
#' @keywords internal
#' @noRd
fgee_qreml_prepare <- function(working, fit.initial, exact = FALSE,
                                rank_tol = 1e-9) {
  if (!inherits(working, "fgee_working_stats")) {
    stop("working must be an fgee_working_stats object.")
  }
  if (working$N < 2L) stop("At least two clusters are required.")

  beta0 <- as.numeric(fit.initial$coefficients)
  H <- 0.5 * (working$W_bar + t(working$W_bar))
  dbar <- as.numeric(working$d_bar)
  b <- as.numeric(H %*% beta0 + dbar)

  ps <- penalty_setup(fit.initial, unpenalized = fit.initial$nsdf)
  lambda_frem <- .fgee_initial_lambda(fit.initial, ps)
  q <- length(lambda_frem)
  components_model <- penalty_components_from_setup(ps, q)$S
  components_model <- lapply(components_model, function(S) 0.5 * (S + t(S)))

  # The ordinary one-step update uses H + P(lambda).  Exact Gaussian fitting
  # uses sum-scale penalized GLS, which is equivalent on the average-cluster
  # scale to H + P(lambda) / N.  Use the corresponding effective penalty in
  # the working criterion while returning the package-scale penalty matrix.
  penalty_scale <- if (isTRUE(exact)) 1 / working$N else 1
  components <- lapply(components_model, `*`, penalty_scale)

  Ssum <- Reduce(`+`, components)
  ee <- eigen(0.5 * (Ssum + t(Ssum)), symmetric = TRUE)
  scale <- max(abs(ee$values), 1)
  keep <- ee$values > rank_tol * scale
  if (!any(keep)) stop("The smoothing penalty has zero numerical rank.")

  Upen <- ee$vectors[, keep, drop = FALSE]
  Sr <- lapply(components, function(S) {
    out <- crossprod(Upen, S %*% Upen)
    0.5 * (out + t(out))
  })

  list(
    working = working,
    fit.initial = fit.initial,
    N = working$N,
    p = length(beta0),
    q = q,
    beta0 = beta0,
    H = H,
    dbar = dbar,
    b = b,
    ps = ps,
    S = components,
    S_model = components_model,
    penalty_scale = penalty_scale,
    exact = isTRUE(exact),
    Upen = Upen,
    Sr = Sr,
    rank_pen = ncol(Upen),
    rank_tol = rank_tol,
    lambda_frem = lambda_frem
  )
}

#' Estimate working-information mismatch
#' @keywords internal
#' @noRd
fgee_qreml_information <- function(prep, center_scores = TRUE) {
  J <- .fgee_working_J(prep$working, center = center_scores)
  H <- prep$H
  fac <- .fgee_factor_pd(H)
  Rinv <- backsolve(fac$R, diag(prep$p))
  M <- crossprod(Rinv, J %*% Rinv)
  M <- 0.5 * (M + t(M))

  phi_all <- mean(diag(M))
  Ssum <- Reduce(`+`, prep$S)
  Sw <- crossprod(Rinv, Ssum %*% Rinv)
  Sw <- 0.5 * (Sw + t(Sw))
  es <- eigen(Sw, symmetric = TRUE)
  keep <- es$values > prep$rank_tol * max(abs(es$values), 1)

  if (any(keep)) {
    Qp <- es$vectors[, keep, drop = FALSE]
    Mpen <- crossprod(Qp, M %*% Qp)
    Mpen <- 0.5 * (Mpen + t(Mpen))
    eig_pen <- pmax(eigen(Mpen, symmetric = TRUE, only.values = TRUE)$values, 0)
    phi_pen <- mean(diag(Mpen))
    anisotropy <- if (is.finite(phi_pen) && phi_pen > 0) {
      norm(Mpen / phi_pen - diag(nrow(Mpen)), type = "F") / sqrt(nrow(Mpen))
    } else NA_real_
    effective_rank <- if (sum(eig_pen^2) > 0) sum(eig_pen)^2 / sum(eig_pen^2) else 0
  } else {
    Mpen <- matrix(0, 0L, 0L)
    eig_pen <- numeric(0)
    phi_pen <- anisotropy <- effective_rank <- NA_real_
  }

  phi_component <- vapply(prep$S, function(Sj) {
    Sjw <- crossprod(Rinv, Sj %*% Rinv)
    den <- sum(diag(Sjw))
    if (!is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
    sum(diag(Sjw %*% M)) / den
  }, numeric(1))

  qfun <- function(x, p) {
    if (!length(x)) return(NA_real_)
    unname(stats::quantile(x, p, names = FALSE, type = 8))
  }

  list(
    J = J,
    H = H,
    M = M,
    M_pen = Mpen,
    phi_all = as.numeric(phi_all),
    phi_pen = as.numeric(phi_pen),
    phi_component = phi_component,
    eig_pen = eig_pen,
    eig_pen_q10 = qfun(eig_pen, 0.10),
    eig_pen_median = if (length(eig_pen)) stats::median(eig_pen) else NA_real_,
    eig_pen_q90 = qfun(eig_pen, 0.90),
    anisotropy = as.numeric(anisotropy),
    effective_rank = as.numeric(effective_rank),
    penalized_rank = sum(keep),
    H_ridge = fac$ridge,
    center_scores = center_scores
  )
}

#' Evaluate qREML and analytic gradient
#' @keywords internal
#' @noRd
fgee_qreml_score_grad <- function(prep, lambda, n_eff,
                                  need_gradient = TRUE) {
  lambda <- as.numeric(lambda)
  if (length(lambda) != prep$q || any(!is.finite(lambda)) || any(lambda <= 0)) {
    return(list(score = Inf, gradient = rep(NA_real_, prep$q)))
  }

  P <- matrix(0, prep$p, prep$p)
  Pr <- matrix(0, prep$rank_pen, prep$rank_pen)
  for (j in seq_len(prep$q)) {
    P <- P + lambda[j] * prep$S[[j]]
    Pr <- Pr + lambda[j] * prep$Sr[[j]]
  }
  P <- 0.5 * (P + t(P))
  Pr <- 0.5 * (Pr + t(Pr))

  fa <- .fgee_factor_pd(prep$H + P)
  fp <- .fgee_factor_pd(Pr)
  beta <- as.numeric(fa$solve(prep$b))
  quad <- sum(prep$b * beta)
  score <- 0.5 * (fa$logdet - fp$logdet - n_eff * quad)

  if (!need_gradient) {
    return(list(
      score = score,
      beta = beta,
      penalty = P,
      n_eff = n_eff,
      ridge_A = fa$ridge,
      ridge_P = fp$ridge
    ))
  }

  grad <- numeric(prep$q)
  for (j in seq_len(prep$q)) {
    dP <- log(10) * lambda[j] * prep$S[[j]]
    dPr <- log(10) * lambda[j] * prep$Sr[[j]]
    trA <- sum(diag(fa$solve(dP)))
    trP <- sum(diag(fp$solve(dPr)))
    betaPbeta <- sum(beta * as.numeric(dP %*% beta))
    grad[j] <- 0.5 * (trA - trP + n_eff * betaPbeta)
  }

  list(
    score = score,
    gradient = grad,
    beta = beta,
    penalty = P,
    n_eff = n_eff,
    ridge_A = fa$ridge,
    ridge_P = fp$ridge
  )
}

#' Optimize sandwich-scaled working qREML
#' @keywords internal
#' @noRd
fgee_tune_qreml <- function(prep,
                             phi_method = c("penalized", "all", "fixed"),
                             phi_fixed = NULL,
                             phi_weight = 1,
                             phi_clip = c(1, 8),
                             center_scores = TRUE,
                             bound = 6,
                             coarse = c(-4, -2, 0, 2, 4),
                             n_starts = 4L,
                             maxit = 100L,
                             factr = 1e8,
                             pgtol = 1e-7,
                             verbose = FALSE) {
  phi_method <- match.arg(phi_method)
  info <- fgee_qreml_information(prep, center_scores = center_scores)
  phi_raw <- switch(
    phi_method,
    penalized = info$phi_pen,
    all = info$phi_all,
    fixed = as.numeric(phi_fixed)[1L]
  )
  if (!is.finite(phi_raw) || phi_raw <= 0) phi_raw <- 1

  phi_weight <- min(max(as.numeric(phi_weight)[1L], 0), 1)
  phi_used <- 1 + phi_weight * (phi_raw - 1)
  if (!is.null(phi_clip)) {
    phi_clip <- sort(as.numeric(phi_clip))
    if (length(phi_clip) != 2L || any(!is.finite(phi_clip)) || phi_clip[1L] <= 0) {
      stop("phi_clip must be NULL or a positive finite length-two vector.")
    }
    phi_used <- min(max(phi_used, phi_clip[1L]), phi_clip[2L])
  }
  n_eff <- prep$N / phi_used

  base <- prep$lambda_frem
  base[!is.finite(base) | base <= 0] <- 1
  starts <- .fgee_qreml_start_grid(prep$q, coarse, bound)

  eval_u <- function(u, grad = TRUE) {
    fgee_qreml_score_grad(prep, base * 10^as.numeric(u), n_eff,
                          need_gradient = grad)
  }
  start_score <- vapply(seq_len(nrow(starts)), function(i) {
    eval_u(starts[i, ], FALSE)$score
  }, numeric(1))
  chosen <- order(start_score)[seq_len(min(n_starts, length(start_score)))]

  n_eval <- 0L
  opts <- vector("list", length(chosen))
  for (s in seq_along(chosen)) {
    last_u <- last <- NULL
    eval_both <- function(u) {
      uu <- as.numeric(u)
      if (!is.null(last_u) && identical(uu, last_u)) return(last)
      last_u <<- uu
      last <<- eval_u(uu, TRUE)
      n_eval <<- n_eval + 1L
      last
    }
    fn <- function(u) {
      z <- eval_both(u)$score
      if (is.finite(z)) z else .Machine$double.xmax^0.25
    }
    gr <- function(u) {
      z <- eval_both(u)$gradient
      z[!is.finite(z)] <- 0
      z
    }
    opts[[s]] <- stats::optim(
      starts[chosen[s], ], fn = fn, gr = gr, method = "L-BFGS-B",
      lower = rep(-bound, prep$q), upper = rep(bound, prep$q),
      control = list(maxit = as.integer(maxit), factr = factr, pgtol = pgtol)
    )
  }

  best <- opts[[which.min(vapply(opts, `[[`, numeric(1), "value"))]]
  ib <- which.min(start_score)
  if (start_score[ib] < best$value) {
    best$par <- starts[ib, ]
    best$value <- start_score[ib]
    best$convergence <- 0L
    best$message <- "best start retained"
  }

  lambda <- base * 10^as.numeric(best$par)
  final <- fgee_qreml_score_grad(prep, lambda, n_eff, need_gradient = TRUE)
  out <- list(
    method = "sandwich_qreml",
    lambda = as.numeric(lambda),
    lambda.star = matrix(as.numeric(lambda), nrow = 1L),
    score = final$score,
    gradient = final$gradient,
    penalty_mat = penalty_from_setup(prep$ps, lambda),
    penalty_effective = final$penalty,
    beta = final$beta,
    phi_method = phi_method,
    phi_raw = phi_raw,
    phi_used = phi_used,
    n_eff = n_eff,
    information = info,
    evaluations = nrow(starts) + n_eval,
    probe_evaluations = nrow(starts),
    gradient_evaluations = n_eval,
    convergence = best$convergence,
    message = best$message,
    boundary = any(abs(best$par) > bound - 0.05),
    log10_multiplier = as.numeric(best$par),
    optim = best,
    prep = prep
  )
  if (verbose) {
    message(
      "sandwich qREML score=", signif(out$score, 8),
      "; phi=", signif(phi_used, 5),
      "; n_eff=", signif(n_eff, 5),
      "; evaluations=", out$evaluations
    )
  }
  out
}
