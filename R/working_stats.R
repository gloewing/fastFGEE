# Compact working-statistics engine -------------------------------------------
#
# This file provides the optimized one-pass replacement for separate W and D
# construction.  The public package still exposes W.estimate() and D.estimate()
# as compatibility wrappers; fgee() uses fgee_build_working_stats() when the
# optimized engine is selected.

.fgee_normalize_corr <- function(x, allow_fpca = FALSE) {
  x <- tolower(as.character(x)[1L])
  if (identical(x, "independence")) x <- "independent"
  allowed <- c("independent", "exchangeable", "ar1")
  if (isTRUE(allow_fpca)) allowed <- c(allowed, "fpca")
  if (!x %in% allowed) {
    stop("Unsupported correlation '", x, "'. Allowed: ",
         paste(allowed, collapse = ", "), call. = FALSE)
  }
  x
}

.fgee_as_rhs_matrix <- function(x) {
  was_vector <- is.null(dim(x))
  if (was_vector) x <- matrix(as.numeric(x), ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  list(x = x, was_vector = was_vector)
}

# Guarded constructor for the SuperGauss Toeplitz backend --------------------
#
# SuperGauss is suggested rather than required: it is reached only through
# corr.solver = "supergauss", or as the fallback when a direct exact operator
# fails, and it is the only dependency that needs a system library (FFTW).
# A default fit with regular-grid AR(1) or exchangeable correlation never calls
# it.  This wrapper keeps the dependency optional while failing informatively
# when the backend is genuinely requested.
#
#' @keywords internal
#' @noRd
.fgee_supergauss_toeplitz <- function(N, acf) {
  if (!requireNamespace("SuperGauss", quietly = TRUE)) {
    stop(
      "The SuperGauss backend requires the 'SuperGauss' package, which is ",
      "suggested rather than required. Install it with ",
      "install.packages(\"SuperGauss\"), or use corr.solver = \"exact\", ",
      "which covers regular-grid AR(1) and exchangeable correlation directly.",
      call. = FALSE
    )
  }
  SuperGauss::Toeplitz$new(N = N, acf = acf)
}

.fgee_apply_exchangeable_inverse <- function(rhs, rho, tol = 1e-10) {
  z <- .fgee_as_rhs_matrix(rhs)
  B <- z$x
  n <- nrow(B)
  if (n <= 1L) return(if (z$was_vector) as.vector(B) else B)

  rho <- as.numeric(rho)[1L]
  d1 <- 1 - rho
  d2 <- 1 + (n - 1) * rho
  if (!is.finite(rho) || d1 <= tol || d2 <= tol) {
    stop("Exchangeable correlation is not positive definite for n=", n,
         " and rho=", format(rho), ".", call. = FALSE)
  }

  cs <- colSums(B)
  out <- B / d1 - (rho / (d1 * d2)) *
    matrix(rep(cs, each = n), nrow = n, ncol = ncol(B))
  if (z$was_vector) as.vector(out) else out
}

.fgee_apply_ar1_regular_inverse <- function(rhs, rho, tol = 1e-10) {
  z <- .fgee_as_rhs_matrix(rhs)
  B <- z$x
  n <- nrow(B)
  if (n <= 1L) return(if (z$was_vector) as.vector(B) else B)

  rho <- as.numeric(rho)[1L]
  den <- 1 - rho^2
  if (!is.finite(rho) || den <= tol) {
    stop("Regular AR(1) requires |rho| < 1; got rho=", format(rho), ".",
         call. = FALSE)
  }

  out <- matrix(0, nrow = n, ncol = ncol(B))
  out[1L, ] <- B[1L, ] - rho * B[2L, ]
  out[n, ] <- B[n, ] - rho * B[n - 1L, ]
  if (n > 2L) {
    ii <- 2L:(n - 1L)
    out[ii, ] <- (1 + rho^2) * B[ii, , drop = FALSE] -
      rho * (B[ii - 1L, , drop = FALSE] + B[ii + 1L, , drop = FALSE])
  }
  out <- out / den
  if (z$was_vector) as.vector(out) else out
}

.fgee_is_regular_grid <- function(times, grid_type = c("auto", "regular", "irregular"),
                                  tol = 1e-10) {
  grid_type <- match.arg(grid_type)
  if (grid_type == "regular") return(TRUE)
  if (grid_type == "irregular") return(FALSE)
  times <- as.numeric(times)
  if (length(times) <= 2L) return(TRUE)
  d <- diff(times)
  all(is.finite(d)) && (max(d) - min(d) <= tol)
}

.fgee_validate_corr_parameter <- function(corr, rho, n, tol = 1e-10) {
  if (corr == "independent" || n <= 1L) return(invisible(TRUE))
  rho <- as.numeric(rho)[1L]
  if (!is.finite(rho)) stop("rho must be a finite scalar.", call. = FALSE)

  if (corr == "ar1" && 1 - rho^2 <= tol) {
    stop("AR(1) requires |rho| < 1; got rho=", format(rho), ".",
         call. = FALSE)
  }
  if (corr == "exchangeable") {
    if (1 - rho <= tol || 1 + (n - 1L) * rho <= tol) {
      stop(
        "Exchangeable correlation is not positive definite for n=", n,
        " and rho=", format(rho), ".", call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.fgee_supergauss_acf <- function(corr, rho, n) {
  switch(
    corr,
    ar1 = rho^(0:(n - 1L)),
    exchangeable = c(1, rep.int(rho, n - 1L)),
    stop("No default SuperGauss ACF for correlation '", corr, "'.",
         call. = FALSE)
  )
}

.fgee_apply_supergauss_acf_inverse <- function(rhs, acf, algo = "gschur",
                                                tol = 1e-8, cache = NULL,
                                                cache_key = NULL) {
  z <- .fgee_as_rhs_matrix(rhs)
  B <- z$x
  n <- nrow(B)
  if (n <= 1L) return(if (z$was_vector) as.vector(B) else B)

  acf <- as.numeric(acf)
  if (length(acf) != n || any(!is.finite(acf)) || abs(acf[1L] - 1) > tol) {
    stop("acf must be a finite length-n vector beginning at 1.", call. = FALSE)
  }
  if (is.null(cache_key)) {
    cache_key <- paste("sg_acf", n, paste(format(acf, digits = 17), collapse = ","),
                       algo, sep = "|")
  }

  Toep <- NULL
  if (!is.null(cache) && exists(cache_key, envir = cache, inherits = FALSE)) {
    Toep <- get(cache_key, envir = cache, inherits = FALSE)
  }
  if (is.null(Toep)) {
    Toep <- .fgee_supergauss_toeplitz(N = n, acf = acf)
    if (!is.null(cache)) assign(cache_key, Toep, envir = cache)
  }
  out <- Toep$solve(B, method = algo, tol = tol)
  if (z$was_vector) as.vector(out) else as.matrix(out)
}

.fgee_apply_supergauss_inverse <- function(rhs, corr, rho, algo = "gschur",
                                            tol = 1e-8, cache = NULL) {
  n <- if (is.null(dim(rhs))) length(rhs) else nrow(rhs)
  acf <- .fgee_supergauss_acf(corr, rho, n)
  key <- paste("sg", corr, n, format(rho, digits = 17), algo, sep = "|")
  .fgee_apply_supergauss_acf_inverse(
    rhs, acf = acf, algo = algo, tol = tol, cache = cache, cache_key = key
  )
}

.fgee_apply_corr_inverse <- function(rhs,
                                      corr,
                                      rho = NULL,
                                      times = NULL,
                                      grid_type = c("auto", "regular", "irregular"),
                                      solver = c("auto", "exact", "supergauss"),
                                      algo = "gschur",
                                      tol = 1e-8,
                                      cache = NULL) {
  corr <- .fgee_normalize_corr(corr, allow_fpca = FALSE)
  solver <- match.arg(solver)
  grid_type <- match.arg(grid_type)

  if (corr == "independent") return(rhs)
  n <- if (is.null(dim(rhs))) length(rhs) else nrow(rhs)
  if (n <= 1L) return(rhs)

  rho <- as.numeric(rho)[1L]
  .fgee_validate_corr_parameter(corr, rho, n, tol = tol)
  if (is.null(times)) times <- seq_len(n)

  regular <- .fgee_is_regular_grid(times, grid_type = grid_type)

  exact_call <- function() {
    if (corr == "exchangeable") {
      return(.fgee_apply_exchangeable_inverse(rhs, rho = rho, tol = tol))
    }
    if (corr == "ar1" && regular) {
      return(.fgee_apply_ar1_regular_inverse(rhs, rho = rho, tol = tol))
    }
    if (corr == "ar1" && !regular) {
      return(fgee_iar1_apply_precision(rhs, time = times, rho = rho))
    }
    stop("No exact inverse operator is available for correlation '", corr, "'.",
         call. = FALSE)
  }

  if (solver == "exact") return(exact_call())
  if (solver == "supergauss") {
    if (!regular && corr == "ar1") return(exact_call())
    return(.fgee_apply_supergauss_inverse(
      rhs, corr = corr, rho = rho, algo = algo, tol = tol, cache = cache
    ))
  }

  out <- try(exact_call(), silent = TRUE)
  if (!inherits(out, "try-error") && all(is.finite(out))) return(out)
  if (!regular && corr == "ar1") {
    cond <- attr(out, "condition")
    msg <- if (inherits(cond, "condition")) conditionMessage(cond) else as.character(out)
    stop(msg, call. = FALSE)
  }

  .fgee_apply_supergauss_inverse(
    rhs, corr = corr, rho = rho, algo = algo, tol = tol, cache = cache
  )
}

.fgee_make_fpca_inverse_operator <- function(fpca_fn, tol = 1e-8) {
  if (is.null(fpca_fn)) stop("fpca_fn is required for corr_fn='fpca'.")

  Phi_full <- as.matrix(fpca_fn$efunctions)
  evals <- as.numeric(fpca_fn$evalues)
  sigma2 <- max(as.numeric(fpca_fn$sigma2), tol)
  argvals <- fpca_fn$argvals

  keep <- which(is.finite(evals) & evals > 0)
  if (!length(keep)) stop("FPCA eigenvalues are not positive and finite.")
  evals <- evals[keep]
  Phi_full <- Phi_full[, keep, drop = FALSE]
  Lambda_inv <- diag(1 / evals, nrow = length(evals))
  cache <- new.env(parent = emptyenv())

  get_parts <- function(idx_vals) {
    pos <- match(idx_vals, argvals)
    if (anyNA(pos)) stop("Some functional indices are absent from fpca_fn$argvals.")
    key <- paste(pos, collapse = ",")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    Phi <- Phi_full[pos, , drop = FALSE]
    M <- Lambda_inv + crossprod(Phi) / sigma2
    M <- 0.5 * (M + t(M))
    fac <- .fgee_factor_pd(M)
    ans <- list(Phi = Phi, solve_mid = fac$solve)
    assign(key, ans, envir = cache)
    ans
  }

  function(rhs, idx_vals) {
    z <- .fgee_as_rhs_matrix(rhs)
    parts <- get_parts(idx_vals)
    out <- z$x / sigma2 -
      parts$Phi %*% parts$solve_mid(crossprod(parts$Phi, z$x)) / sigma2^2
    if (z$was_vector) as.vector(out) else out
  }
}

.fgee_apply_kron_inverse <- function(rhs, n_fun, n_long,
                                      apply_fun, apply_long) {
  z <- .fgee_as_rhs_matrix(rhs)
  B <- z$x
  q <- ncol(B)
  if (nrow(B) != n_fun * n_long) {
    stop("Kronecker RHS has incompatible dimensions.", call. = FALSE)
  }

  A <- array(B, dim = c(n_fun, n_long, q))
  F <- matrix(A, nrow = n_fun, ncol = n_long * q)
  F <- apply_fun(F)
  A <- array(F, dim = c(n_fun, n_long, q))

  L <- aperm(A, c(2, 1, 3))
  L <- matrix(L, nrow = n_long, ncol = n_fun * q)
  L <- apply_long(L)
  A <- aperm(array(L, dim = c(n_long, n_fun, q)), c(2, 1, 3))

  out <- matrix(A, nrow = n_fun * n_long, ncol = q)
  if (z$was_vector) as.vector(out) else out
}

.fgee_unique_scalar <- function(x, label) {
  u <- unique(as.numeric(x))
  u <- u[is.finite(u)]
  if (length(u) != 1L) stop(label, " must be a finite scalar within this block.")
  u
}


.fgee_tensor_grid_info <- function(idx_fn, idx_long, index_fn, index_long,
                                   require_complete = FALSE) {
  n <- length(idx_fn)
  if (length(idx_long) != n) stop("Functional and longitudinal indices differ in length.")

  fn_levels <- sort(unique(idx_fn))
  long_levels <- sort(unique(idx_long))
  n_fun <- length(fn_levels)
  n_long <- length(long_levels)

  fn_pos <- match(idx_fn, fn_levels)
  long_pos <- match(idx_long, long_levels)
  key <- fn_pos + (long_pos - 1L) * n_fun
  duplicated_cell <- anyDuplicated(key) > 0L
  complete <- !duplicated_cell && n == n_fun * n_long &&
    length(unique(key)) == n_fun * n_long

  if (isTRUE(require_complete) && !complete) {
    if (duplicated_cell) {
      stop(
        "Kronecker path requires exactly one row per (", index_long, ", ",
        index_fn, ") cell within each cluster.", call. = FALSE
      )
    }
    stop(
      "Kronecker path requires a complete tensor-product grid within each cluster: ",
      "all combinations of (", index_long, ", ", index_fn, ") must be present.",
      call. = FALSE
    )
  }

  list(
    complete = complete,
    duplicated = duplicated_cell,
    n_fun = n_fun,
    n_long = n_long,
    fn_levels = fn_levels,
    long_levels = long_levels
  )
}

.fgee_apply_cluster_inverse <- function(rhs,
                                         idx_fn,
                                         idx_long,
                                         corr_fn,
                                         corr_long,
                                         rho_fn = NULL,
                                         rho_long = NULL,
                                         fpca_apply = NULL,
                                         grid_type = "auto",
                                         corr_solver = "auto",
                                         algo = "gschur",
                                         tol = 1e-8,
                                         cache = NULL,
                                         index_fn = "yindex.vec",
                                         index_long = "time") {
  corr_fn <- .fgee_normalize_corr(corr_fn, allow_fpca = TRUE)
  corr_long <- .fgee_normalize_corr(corr_long, allow_fpca = FALSE)
  z <- .fgee_as_rhs_matrix(rhs)
  Q <- z$x

  if (corr_fn == "independent" && corr_long == "independent") {
    return(if (z$was_vector) as.vector(Q) else Q)
  }

  both <- corr_fn != "independent" && corr_long != "independent"
  grid <- .fgee_tensor_grid_info(
    idx_fn, idx_long, index_fn = index_fn, index_long = index_long,
    require_complete = both
  )

  if (both) {
    n_fun <- grid$n_fun
    n_long <- grid$n_long
    times_fn <- grid$fn_levels
    times_long <- grid$long_levels

    if (corr_fn == "fpca") {
      if (is.null(fpca_apply)) stop("FPCA inverse operator is unavailable.")
      af <- function(B) fpca_apply(B, times_fn)
    } else {
      rf <- .fgee_unique_scalar(rho_fn, "rho_fn")
      af <- function(B) .fgee_apply_corr_inverse(
        B, corr = corr_fn, rho = rf, times = times_fn,
        grid_type = grid_type, solver = corr_solver,
        algo = algo, tol = tol, cache = cache
      )
    }

    rl <- .fgee_unique_scalar(rho_long, "rho_long")
    al <- function(B) .fgee_apply_corr_inverse(
      B, corr = corr_long, rho = rl, times = times_long,
      grid_type = grid_type, solver = corr_solver,
      algo = algo, tol = tol, cache = cache
    )

    if (corr_fn != "fpca") {
      k_out <- .fgee_kron_kernel_try(
        Q, n_fun, n_long, corr_fn, rf, corr_long, rl,
        times_fn = times_fn, times_long = times_long,
        grid_type = grid_type, solver = corr_solver, tol = tol
      )
      if (!is.null(k_out)) {
        return(if (z$was_vector) as.vector(k_out) else k_out)
      }
    }
    out <- .fgee_apply_kron_inverse(Q, n_fun, n_long, af, al)
    return(if (z$was_vector) as.vector(out) else out)
  }

  out <- matrix(NA_real_, nrow = nrow(Q), ncol = ncol(Q))

  if (corr_long != "independent") {
    rho_unique <- unique(as.numeric(rho_long[is.finite(rho_long)]))
    if (isTRUE(grid$complete) && length(rho_unique) == 1L) {
      k_out <- .fgee_kron_kernel_try(
        Q, grid$n_fun, grid$n_long, "independent", 0, corr_long, rho_unique,
        times_fn = grid$fn_levels, times_long = grid$long_levels,
        grid_type = grid_type, solver = corr_solver, tol = tol
      )
      if (!is.null(k_out)) {
        return(if (z$was_vector) as.vector(k_out) else k_out)
      }
      A <- array(Q, dim = c(grid$n_fun, grid$n_long, ncol(Q)))
      L <- aperm(A, c(2, 1, 3))
      L <- matrix(L, nrow = grid$n_long,
                  ncol = grid$n_fun * ncol(Q))
      L <- .fgee_apply_corr_inverse(
        L, corr = corr_long, rho = rho_unique,
        times = grid$long_levels, grid_type = grid_type,
        solver = corr_solver, algo = algo, tol = tol, cache = cache
      )
      A <- aperm(
        array(L, dim = c(grid$n_long, grid$n_fun, ncol(Q))),
        c(2, 1, 3)
      )
      out <- matrix(A, nrow = nrow(Q), ncol = ncol(Q))
    } else {
      groups <- split(seq_along(idx_fn), idx_fn, drop = TRUE)
      for (ii in groups) {
        rr <- .fgee_unique_scalar(rho_long[ii], "rho_long")
        out[ii, ] <- .fgee_apply_corr_inverse(
          Q[ii, , drop = FALSE], corr = corr_long, rho = rr,
          times = idx_long[ii], grid_type = grid_type,
          solver = corr_solver, algo = algo, tol = tol, cache = cache
        )
      }
    }
  } else {
    rho_unique <- if (corr_fn == "fpca") numeric(0) else {
      unique(as.numeric(rho_fn[is.finite(rho_fn)]))
    }
    can_batch <- isTRUE(grid$complete) &&
      (corr_fn == "fpca" || length(rho_unique) == 1L)

    if (can_batch) {
      if (corr_fn != "fpca") {
        k_out <- .fgee_kron_kernel_try(
          Q, grid$n_fun, grid$n_long, corr_fn, rho_unique, "independent", 0,
          times_fn = grid$fn_levels, times_long = grid$long_levels,
          grid_type = grid_type, solver = corr_solver, tol = tol
        )
        if (!is.null(k_out)) {
          return(if (z$was_vector) as.vector(k_out) else k_out)
        }
      }
      A <- array(Q, dim = c(grid$n_fun, grid$n_long, ncol(Q)))
      F <- matrix(A, nrow = grid$n_fun,
                  ncol = grid$n_long * ncol(Q))
      if (corr_fn == "fpca") {
        if (is.null(fpca_apply)) stop("FPCA inverse operator is unavailable.")
        F <- fpca_apply(F, grid$fn_levels)
      } else {
        F <- .fgee_apply_corr_inverse(
          F, corr = corr_fn, rho = rho_unique,
          times = grid$fn_levels, grid_type = grid_type,
          solver = corr_solver, algo = algo, tol = tol, cache = cache
        )
      }
      out <- matrix(
        array(F, dim = c(grid$n_fun, grid$n_long, ncol(Q))),
        nrow = nrow(Q), ncol = ncol(Q)
      )
    } else {
      groups <- split(seq_along(idx_long), idx_long, drop = TRUE)
      for (ii in groups) {
        if (corr_fn == "fpca") {
          if (is.null(fpca_apply)) stop("FPCA inverse operator is unavailable.")
          out[ii, ] <- fpca_apply(Q[ii, , drop = FALSE], idx_fn[ii])
        } else {
          rr <- .fgee_unique_scalar(rho_fn[ii], "rho_fn")
          out[ii, ] <- .fgee_apply_corr_inverse(
            Q[ii, , drop = FALSE], corr = corr_fn, rho = rr,
            times = idx_fn[ii], grid_type = grid_type,
            solver = corr_solver, algo = algo, tol = tol, cache = cache
          )
        }
      }
    }
  }

  if (z$was_vector) as.vector(out) else out
}

.fgee_match_cluster_fold <- function(ids, fold_id) {
  if (is.null(fold_id)) return(NULL)
  if (!is.null(names(fold_id))) {
    if (!all(as.character(ids) %in% names(fold_id))) {
      stop("Named fold_id does not include every cluster.")
    }
    fold_id <- fold_id[as.character(ids)]
  }
  fold_id <- as.integer(fold_id)
  if (length(fold_id) != length(ids) || anyNA(fold_id) || any(fold_id < 1L)) {
    stop("fold_id must contain one positive integer per cluster.")
  }
  fold_id
}

#' Build compact cluster-level working statistics
#' @keywords internal
#' @noRd
.fgee_build_working_stats_core <- function(
    dx,
    namesd,
    cname_ = "cname_",
    corr_fn = "independent",
    corr_long = "independent",
    index_fn = "yindex.vec",
    index_long = "time",
    resid_col = "resid",
    fpca_fn = NULL,
    retain = c("scores", "aggregate", "full"),
    corr_solver = c("auto", "exact", "supergauss"),
    algo = "gschur",
    grid_type = c("auto", "regular", "irregular"),
    tol = 1e-8,
    beta0 = NULL,
    exact_gaussian = FALSE,
    gaussian_fold_id = NULL,
    copy_dt = TRUE,
    ensure_order = TRUE,
    check = TRUE) {

  retain <- match.arg(retain)
  corr_solver <- match.arg(corr_solver)
  grid_type <- match.arg(grid_type)
  corr_fn <- .fgee_normalize_corr(corr_fn, allow_fpca = TRUE)
  corr_long <- .fgee_normalize_corr(corr_long, allow_fpca = FALSE)

  dd <- data.table::as.data.table(dx)
  if (isTRUE(copy_dt)) dd <- data.table::copy(dd)

  required <- unique(c(
    namesd, cname_, index_fn, index_long, resid_col, "muprime",
    if (isTRUE(exact_gaussian) || !is.null(gaussian_fold_id)) "Y" else character(0)
  ))
  missing <- setdiff(required, names(dd))
  if (length(missing)) stop("Missing working-data columns: ", paste(missing, collapse = ", "))
  if (!"sqrtv" %in% names(dd)) {
    if (!"v" %in% names(dd)) stop("Need column 'v' or 'sqrtv'.")
    dd[, sqrtv := sqrt(v)]
  }
  if (corr_fn %in% c("ar1", "exchangeable") && !"rho_fn" %in% names(dd)) {
    stop("Need rho_fn column for corr_fn='", corr_fn, "'.")
  }
  if (corr_long != "independent" && !"rho_long" %in% names(dd)) {
    stop("Need rho_long column for corr_long='", corr_long, "'.")
  }

  if (isTRUE(ensure_order)) {
    data.table::setorderv(dd, c(cname_, index_long, index_fn))
  }

  cid <- as.character(dd[[cname_]])
  rr <- rle(cid)
  ids <- rr$values
  if (anyDuplicated(ids)) {
    stop(
      "Cluster rows are not contiguous. Sort by '", cname_, "', '",
      index_long, "', and '", index_fn, "' or set ensure_order=TRUE.",
      call. = FALSE
    )
  }
  ends <- cumsum(rr$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  if (isTRUE(check) && !isTRUE(ensure_order)) {
    bad <- vapply(seq_along(starts), function(i) {
      ii <- starts[i]:ends[i]
      il <- dd[[index_long]][ii]
      jf <- dd[[index_fn]][ii]
      if (length(ii) <= 1L) return(FALSE)
      pil <- head(il, -1L)
      pjf <- head(jf, -1L)
      any(tail(il, -1L) < pil |
            (tail(il, -1L) == pil & tail(jf, -1L) < pjf))
    }, logical(1))
    if (any(bad)) {
      stop(
        "Rows are not in canonical (cluster, longitudinal, functional) order; ",
        "set ensure_order=TRUE.", call. = FALSE
      )
    }
  }
  N <- length(ids)
  p <- length(namesd)

  fold_id <- .fgee_match_cluster_fold(ids, gaussian_fold_id)
  gaussian_cv <- NULL
  if (!is.null(fold_id)) {
    K <- max(fold_id)
    gaussian_cv <- list(
      C = array(0, dim = c(p, p, K)),
      c = matrix(0, nrow = p, ncol = K),
      s = numeric(K),
      fold_id = setNames(fold_id, ids),
      fold_cluster_count = tabulate(fold_id, nbins = K)
    )
  }

  D <- if (retain != "aggregate") matrix(NA_real_, p, N) else NULL
  D_exact <- if (isTRUE(exact_gaussian) && retain != "aggregate") {
    matrix(NA_real_, p, N)
  } else NULL
  W_list <- if (retain == "full") vector("list", N) else NULL

  W_sum <- matrix(0, p, p)
  d_sum <- numeric(p)
  dd_sum <- matrix(0, p, p)
  d_exact_sum <- if (isTRUE(exact_gaussian)) numeric(p) else NULL
  cluster_size <- rr$lengths

  corr_cache <- new.env(parent = emptyenv())
  fpca_apply <- if (corr_fn == "fpca") {
    .fgee_make_fpca_inverse_operator(fpca_fn, tol = tol)
  } else NULL

  for (i in seq_len(N)) {
    ii <- starts[i]:ends[i]
    X <- as.matrix(dd[ii, ..namesd])
    storage.mode(X) <- "double"
    sv <- as.numeric(dd$sqrtv[ii])
    mp <- as.numeric(dd$muprime[ii])
    resid <- as.numeric(dd[[resid_col]][ii])

    if (isTRUE(check)) {
      if (any(!is.finite(X)) || any(!is.finite(sv)) || any(sv <= 0) ||
          any(!is.finite(mp)) || any(!is.finite(resid))) {
        stop("Non-finite working quantities in cluster '", ids[i], "'.")
      }
    }

    Z <- X * (mp / sv)
    Q <- cbind(Z, resid)
    rho_fn_i <- if ("rho_fn" %in% names(dd)) dd$rho_fn[ii] else NULL
    rho_long_i <- if ("rho_long" %in% names(dd)) dd$rho_long[ii] else NULL

    # Wi and di are the only consumers of R^{-1} Q, so when both axes are
    # exchangeable or independent the Gram matrix can be formed directly from
    # rank-one sufficient statistics without applying the inverse at all.
    # NULL falls through to the reference apply path below.
    gram <- .fgee_cluster_gram(
      Q,
      idx_fn = dd[[index_fn]][ii],
      idx_long = dd[[index_long]][ii],
      corr_fn = corr_fn,
      corr_long = corr_long,
      rho_fn = rho_fn_i,
      rho_long = rho_long_i,
      grid_type = grid_type,
      corr_solver = corr_solver,
      index_fn = index_fn,
      index_long = index_long
    )

    if (is.null(gram)) {
      rinv_Q <- .fgee_apply_cluster_inverse(
        Q,
        idx_fn = dd[[index_fn]][ii],
        idx_long = dd[[index_long]][ii],
        corr_fn = corr_fn,
        corr_long = corr_long,
        rho_fn = rho_fn_i,
        rho_long = rho_long_i,
        fpca_apply = fpca_apply,
        grid_type = grid_type,
        corr_solver = corr_solver,
        algo = algo,
        tol = tol,
        cache = corr_cache,
        index_fn = index_fn,
        index_long = index_long
      )
      Wi <- crossprod(Z, rinv_Q[, seq_len(p), drop = FALSE])
      di <- as.numeric(crossprod(Z, rinv_Q[, p + 1L]))
    } else {
      Wi <- gram[seq_len(p), seq_len(p), drop = FALSE]
      di <- as.numeric(gram[seq_len(p), p + 1L])
    }
    Wi <- 0.5 * (Wi + t(Wi))

    if (any(!is.finite(Wi)) || any(!is.finite(di))) {
      stop("Non-finite W or d in cluster '", ids[i], "'.")
    }

    W_sum <- W_sum + Wi
    d_sum <- d_sum + di
    dd_sum <- dd_sum + tcrossprod(di)

    if (!is.null(D)) D[, i] <- di
    if (!is.null(W_list)) W_list[[i]] <- Wi

    if (isTRUE(exact_gaussian)) {
      if (is.null(beta0) || length(beta0) != p) {
        stop("exact_gaussian=TRUE requires beta0 with length p.")
      }
      dex <- di + as.numeric(Wi %*% beta0)
      d_exact_sum <- d_exact_sum + dex
      if (!is.null(D_exact)) D_exact[, i] <- dex
    }

    if (!is.null(gaussian_cv)) {
      k <- fold_id[i]
      mi <- length(ii)
      sc <- 1 / (N * mi)
      yi <- as.numeric(dd$Y[ii])
      gaussian_cv$C[, , k] <- gaussian_cv$C[, , k] + sc * crossprod(X)
      gaussian_cv$c[, k] <- gaussian_cv$c[, k] + sc * as.numeric(crossprod(X, yi))
      gaussian_cv$s[k] <- gaussian_cv$s[k] + sc * sum(yi^2)
    }
  }

  if (!is.null(D)) {
    dimnames(D) <- list(namesd, ids)
    D <- unname(D)
    colnames(D) <- ids
    rownames(D) <- namesd
  }
  if (!is.null(D_exact)) {
    colnames(D_exact) <- ids
    rownames(D_exact) <- namesd
  }
  if (!is.null(W_list)) names(W_list) <- ids

  out <- list(
    cluster_id = ids,
    cluster_size = as.integer(cluster_size),
    N = N,
    p = p,
    coefficient_names = namesd,
    W_sum = W_sum,
    W_bar = W_sum / N,
    d_sum = d_sum,
    d_bar = d_sum / N,
    dd_sum = dd_sum,
    D = D,
    W = W_list,
    D_exact = D_exact,
    d_exact_sum = d_exact_sum,
    gaussian_cv = gaussian_cv,
    retain = retain,
    corr_fn = corr_fn,
    corr_long = corr_long,
    corr_solver = corr_solver,
    algo = algo,
    grid_type = grid_type,
    cname_ = cname_,
    index_fn = index_fn,
    index_long = index_long
  )
  class(out) <- c("fgee_working_stats", "list")
  out
}

#' @keywords internal
#' @noRd
print.fgee_working_stats <- function(x, ...) {
  cat("<fgee_working_stats>\n")
  cat("  clusters:", x$N, "\n")
  cat("  coefficients:", x$p, "\n")
  cat("  retention:", x$retain, "\n")
  cat("  correlation:", x$corr_long, "x", x$corr_fn, "\n")
  invisible(x)
}

.fgee_working_dmat <- function(x, exact = FALSE, required = TRUE) {
  if (!inherits(x, "fgee_working_stats")) stop("Expected fgee_working_stats.")
  ans <- if (isTRUE(exact)) x$D_exact else x$D
  if (is.null(ans) && isTRUE(required)) {
    stop("Cluster scores were not retained. Use working.retain='scores' or 'full'.")
  }
  ans
}

.fgee_working_dlist <- function(x, exact = FALSE, required = TRUE) {
  D <- .fgee_working_dmat(x, exact = exact, required = required)
  if (is.null(D)) return(NULL)
  ans <- lapply(seq_len(ncol(D)), function(i) as.numeric(D[, i]))
  names(ans) <- x$cluster_id
  ans
}

.fgee_working_wlist <- function(x, required = TRUE) {
  if (!inherits(x, "fgee_working_stats")) stop("Expected fgee_working_stats.")
  if (is.null(x$W) && isTRUE(required)) {
    stop("Cluster W_i matrices were not retained. Use working.retain='full'.")
  }
  x$W
}

.fgee_working_J <- function(x, center = TRUE) {
  if (!inherits(x, "fgee_working_stats")) stop("Expected fgee_working_stats.")
  if (x$N < 2L) stop("At least two clusters are required.")
  if (isTRUE(center)) {
    J <- (x$dd_sum - x$N * tcrossprod(x$d_bar)) / (x$N - 1L)
  } else {
    J <- x$dd_sum / x$N
  }
  0.5 * (J + t(J))
}

.fgee_resolve_working_retain <- function(requested = c("auto", "scores", "aggregate", "full"),
                                          joint.CI = "wild",
                                          var.type = "sandwich",
                                          exact = FALSE,
                                          max.iter = 1L,
                                          sp.method = "fastk_staged",
                                          tune.method = "one-step",
                                          gee.fit = TRUE) {
  requested <- match.arg(requested)

  need_full <- identical(var.type, "boot") ||
    identical(tune.method, "fully-iterated") ||
    (!identical(joint.CI, FALSE) && !identical(joint.CI, "wild"))

  need_scores <- need_full || identical(joint.CI, "wild") ||
    identical(var.type, "fastboot") ||
    (isTRUE(gee.fit) &&
       sp.method %in% c("fastk_staged", "fastk_grad", "fastk_grad_fast",
                       "qreml_fastk", "auto"))

  needed <- if (need_full) "full" else if (need_scores) "scores" else "aggregate"
  if (requested == "auto") return(needed)

  level <- c(aggregate = 1L, scores = 2L, full = 3L)
  if (level[[requested]] < level[[needed]]) {
    stop(
      "working.retain='", requested, "' is insufficient for this workflow; ",
      "use working.retain='", needed, "' (or 'auto').",
      call. = FALSE
    )
  }
  requested
}

