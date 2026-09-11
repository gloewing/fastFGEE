# Fast start selection for the analytic-gradient fastK selector ---------------
#
# The stock `fgee_tune_fastk_grad()` spends most of its budget before L-BFGS-B
# begins.  For the common q = 3 case it scores a 55-point isotropic ray plus a
# qREML candidate plus 2q axis probes -- about 62 exact evaluations -- and then
# runs four independent optimisations.  Measured evaluation counts of 66 to 108
# therefore decompose as roughly 62 start probes and only 4 to 46 optimiser
# calls, so the start search, not the optimiser, dominates.
#
# This file implements a cheaper initialiser that follows the published staged
# structure at initialisation resolution:
#
#   Stage 1  a coarse common-scale ray, lambda = 10^c * lambda_fREML for
#            c in {-3, ..., 3}: seven evaluations rather than 55.  The grid is
#            extended outward while the best point sits on an endpoint, then
#            optionally refined by one-dimensional Brent search inside the
#            bracketing interval.
#   Stage 2  a sparse anisotropic search around the Stage 1 result, letting the
#            smoothing parameters separate: a 3^q grid for small q, or a
#            coordinate sweep otherwise, so the cost is O(mq) rather than the
#            O(m^q) of a full Cartesian stage.
#   Polish   a single L-BFGS-B run from the best exact-score candidate, with a
#            second run only when diagnostics call for one.
#
# The optimised criterion is unchanged throughout: every candidate is scored
# with the exact full fastK objective and the returned lambda is the argmin over
# everything evaluated.
#
# @keywords internal
# @noRd
NULL

#' Adaptive common-scale (Stage 1) search
#'
#' Locates the order of magnitude of the smoothing parameters by scaling the
#' whole fREML vector by a single multiplier.  Expansion continues while the
#' discrete minimum lies on an endpoint, because stopping early would report the
#' edge of the search window rather than an optimum.
#'
#' @keywords internal
#' @noRd
.fgee_fastk_stage1_adaptive <- function(obj, q, c_grid = -3:3, bound = 8,
                                        max_expand = 32L, brent = TRUE) {
  cs <- sort(unique(as.numeric(c_grid)))
  cs <- cs[abs(cs) <= bound]
  sc <- vapply(cs, function(c) obj$eval(rep(c, q), FALSE)$score, numeric(1))

  for (e in seq_len(max_expand)) {
    i <- which.min(sc)
    step <- if (length(cs) > 1L) min(diff(cs)) else 1
    if (i == 1L && cs[1L] - step >= -bound) {
      newc <- cs[1L] - step
    } else if (i == length(cs) && cs[length(cs)] + step <= bound) {
      newc <- cs[length(cs)] + step
    } else {
      break
    }
    news <- obj$eval(rep(newc, q), FALSE)$score
    cs <- c(cs, newc)
    sc <- c(sc, news)
    o <- order(cs)
    cs <- cs[o]
    sc <- sc[o]
  }

  i <- which.min(sc)
  c_best <- cs[i]
  interior <- (i > 1L && i < length(cs)) || abs(cs[i]) >= bound - 1e-9

  if (isTRUE(brent) && length(cs) >= 3L && i > 1L && i < length(cs)) {
    op <- try(stats::optimize(
      function(c) obj$eval(rep(c, q), FALSE)$score,
      interval = c(cs[i - 1L], cs[i + 1L]), tol = 0.05
    ), silent = TRUE)
    if (!inherits(op, "try-error") && is.finite(op$objective) &&
        op$objective <= sc[i]) {
      c_best <- op$minimum
    }
  }

  list(u = rep(c_best, q), score = obj$eval(rep(c_best, q), FALSE)$score,
       c_best = c_best, grid = cs, grid_score = sc, interior = interior)
}

#' Sparse anisotropic (Stage 2) search
#'
#' Allows the smoothing parameters to separate around the Stage 1 centre without
#' paying for a full Cartesian grid.  A `5^q` stage costs 125 evaluations at
#' q = 3, which is more than an entire gradient run, so a `3^q` grid is used for
#' small q and a coordinate sweep beyond that.
#'
#' @keywords internal
#' @noRd
.fgee_fastk_stage2_sparse <- function(obj, u_center, q, bound = 8,
                                      offsets = c(-1, 0, 1),
                                      max_q_cartesian = 3L,
                                      sweep_offsets = c(-2, -1, 1, 2),
                                      expand_edges = TRUE) {
  clamp <- function(u) pmin(pmax(u, -bound), bound)

  if (q <= max_q_cartesian) {
    grid <- as.matrix(expand.grid(rep(list(offsets), q), KEEP.OUT.ATTRS = FALSE))
    cand <- t(apply(grid, 1L, function(d) clamp(u_center + d)))
    if (q == 1L) cand <- matrix(cand, ncol = 1L)
    sc <- apply(cand, 1L, function(u) obj$eval(u, FALSE)$score)
    i <- which.min(sc)
    return(list(u = cand[i, ], score = sc[i], mode = "cartesian",
                n_points = nrow(cand)))
  }

  u <- clamp(as.numeric(u_center))
  best <- obj$eval(u, FALSE)$score
  n_pts <- 1L
  for (j in seq_len(q)) {
    offs <- sweep_offsets
    repeat {
      cand_j <- clamp(u[j] + offs)
      sc_j <- vapply(cand_j, function(v) {
        uu <- u; uu[j] <- v; obj$eval(uu, FALSE)$score
      }, numeric(1))
      n_pts <- n_pts + length(cand_j)
      m <- which.min(sc_j)
      if (sc_j[m] < best) {
        best <- sc_j[m]
        u[j] <- cand_j[m]
      } else {
        break
      }
      if (!isTRUE(expand_edges) || !(m == 1L || m == length(cand_j))) break
      if (abs(u[j]) >= bound - 1e-9) break
      offs <- offs * 2
    }
  }
  list(u = u, score = best, mode = "coordinate", n_points = n_pts)
}

#' Counting, memoising objective in log10-multiplier space
#'
#' The stock selector re-evaluates each chosen start inside `optim()` after
#' having already scored it during probing, and its one-point cache is reset for
#' every start.  This cache is shared across all stages and all starts, so no
#' point is ever paid for twice.
#'
#' `kernel = TRUE` routes the per-fold loss and score contributions through the
#' compiled kernel when it is available; the result is identical to the R path to
#' within floating-point reassociation (verified to about 1e-16).
#'
#' @keywords internal
#' @noRd
.fgee_fastk_objective <- function(prep, base, kernel = fgee_fastk_kernel_ok(),
                                  digits = 10L) {
  use_k <- isTRUE(kernel) && !is.null(prep$X_eval) && !isTRUE(prep$gaussian) &&
    !is.null(.fgee_fastk_family_code(prep$family_key))
  env <- new.env(parent = emptyenv())
  env$cache <- new.env(parent = emptyenv())
  env$n_value <- 0L
  env$n_grad <- 0L
  key <- function(u) paste(round(as.numeric(u), digits), collapse = ",")

  eval_u <- function(u, need_gradient = FALSE) {
    k <- key(u)
    hit <- env$cache[[k]]
    if (!is.null(hit) && (!need_gradient || !is.null(hit$gradient))) return(hit)
    lambda <- base * 10^as.numeric(u)
    res <- if (use_k) {
      .fgee_fastk_score_grad_kernel(prep, lambda, need_gradient = need_gradient)
    } else {
      fgee_fastk_score_grad(prep, lambda, need_gradient = need_gradient)
    }
    if (need_gradient) env$n_grad <- env$n_grad + 1L else env$n_value <- env$n_value + 1L
    env$cache[[k]] <- res
    res
  }
  list(eval = eval_u,
       counts = function() list(value = env$n_value, gradient = env$n_grad),
       kernel = use_k)
}

#' Fast analytic-gradient exact fastK
#'
#' Adaptive Stage 1 plus sparse Stage 2 plus a single continuous polish.  Returns
#' the same structure as [fgee_tune_fastk_grad()] so the dispatcher and engine
#' consume it unchanged.
#'
#' @keywords internal
#' @noRd
fgee_tune_fastk_grad_fast <- function(prep,
                                      qreml_lambda = NULL,
                                      bound = 8,
                                      c_grid = -3:3,
                                      stage1_brent = TRUE,
                                      stage2_offsets = c(-1, 0, 1),
                                      max_q_cartesian = 3L,
                                      maxit = 100L,
                                      factr = 1e8,
                                      pgtol = 1e-7,
                                      second_start = TRUE,
                                      pg_tol_rel = 1e-3,
                                      kernel = fgee_fastk_kernel_ok(),
                                      verbose = FALSE) {
  base <- prep$lambda_frem
  base[!is.finite(base) | base <= 0] <- 1
  q <- prep$q
  obj <- .fgee_fastk_objective(prep, base, kernel = kernel)
  clamp <- function(u) pmin(pmax(u, -bound), bound)

  cand <- list(frem = rep(0, q))
  if (!is.null(qreml_lambda)) {
    uq <- log10(as.numeric(qreml_lambda) / base)
    if (all(is.finite(uq))) cand$qreml <- clamp(uq)
  }
  s1 <- .fgee_fastk_stage1_adaptive(obj, q, c_grid = c_grid, bound = bound,
                                    brent = stage1_brent)
  cand$stage1 <- s1$u
  s2 <- .fgee_fastk_stage2_sparse(obj, s1$u, q, bound = bound,
                                  offsets = stage2_offsets,
                                  max_q_cartesian = max_q_cartesian)
  cand$stage2 <- s2$u

  cand_mat <- do.call(rbind, cand)
  cand_sc <- apply(cand_mat, 1L, function(u) obj$eval(u, FALSE)$score)
  ord <- order(cand_sc)

  fn <- function(u) {
    v <- obj$eval(u, TRUE)$score
    if (is.finite(v)) v else .Machine$double.xmax^0.25
  }
  gr <- function(u) {
    g <- obj$eval(u, TRUE)$gradient
    g[!is.finite(g)] <- 0
    g
  }
  run_opt <- function(u0) {
    stats::optim(u0, fn = fn, gr = gr, method = "L-BFGS-B",
                 lower = rep(-bound, q), upper = rep(bound, q),
                 control = list(maxit = as.integer(maxit), factr = factr,
                                pgtol = pgtol))
  }

  u_start <- cand_mat[ord[1L], ]
  fit <- run_opt(u_start)

  # Relative projected gradient, scaled by the gradient magnitude at the start.
  # Scaling by the score instead would make this ratio vanishingly small,
  # because the score is O(M), and the trigger below would never fire.
  g_start <- obj$eval(u_start, TRUE)$gradient
  g_scale <- if (all(is.finite(g_start))) max(max(abs(g_start)), 1e-12) else 1
  pg_rel <- function(u) {
    g <- obj$eval(u, TRUE)$gradient
    if (!all(is.finite(g))) return(Inf)
    g[abs(u) >= bound - 1e-8] <- 0
    max(abs(g)) / g_scale
  }

  n_second <- 0L
  if (isTRUE(second_start) && nrow(cand_mat) >= 2L) {
    runner <- cand_mat[ord[2L], ]
    close_score <- is.finite(cand_sc[ord[2L]]) &&
      abs(cand_sc[ord[2L]] - cand_sc[ord[1L]]) <=
        1e-3 * max(abs(cand_sc[ord[1L]]), 1e-12)
    far_lambda <- max(abs(runner - u_start)) > 1
    need_second <- !identical(as.integer(fit$convergence), 0L) ||
      any(abs(fit$par) > bound - 0.05) ||
      pg_rel(fit$par) > pg_tol_rel ||
      (close_score && far_lambda)
    if (isTRUE(need_second)) {
      f2 <- run_opt(runner)
      n_second <- 1L
      if (is.finite(f2$value) && f2$value < fit$value) fit <- f2
    }
  }

  ib <- ord[1L]
  if (cand_sc[ib] < fit$value) {
    fit$par <- cand_mat[ib, ]
    fit$value <- cand_sc[ib]
    fit$convergence <- 0L
    fit$message <- "best candidate retained"
  }

  cnt <- obj$counts()
  lambda <- base * 10^as.numeric(fit$par)
  if (isTRUE(verbose)) {
    message("fast gradient fastK score=", signif(fit$value, 8),
            "; value evaluations=", cnt$value,
            "; gradient evaluations=", cnt$gradient,
            "; kernel=", obj$kernel)
  }

  list(
    method = "fastk_grad_fast",
    lambda = as.numeric(lambda),
    lambda.star = matrix(as.numeric(lambda), nrow = 1L),
    score = as.numeric(fit$value),
    penalty_mat = penalty_from_setup(prep$ps, lambda),
    evaluations = cnt$value + cnt$gradient,
    probe_evaluations = cnt$value,
    gradient_evaluations = cnt$gradient,
    convergence = fit$convergence,
    message = fit$message,
    boundary = any(abs(fit$par) > bound - 0.05),
    log10_multiplier = as.numeric(fit$par),
    optim = fit,
    kernel = obj$kernel,
    stage1_c = s1$c_best,
    stage1_interior = isTRUE(s1$interior),
    stage2_mode = s2$mode,
    stage2_points = s2$n_points,
    chosen_candidate = rownames(cand_mat)[ib],
    n_second_start = n_second,
    workspace = prep
  )
}
