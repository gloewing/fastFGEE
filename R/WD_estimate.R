
# ============================================================================
# .getW()
# ============================================================================
#' @keywords internal
#' @noRd
.getW <- function(dd, namesd, cname_,
                  index_fn = "yindex.vec", index_long = "time",
                  second.deriv = FALSE,
                  corr_fn = "ar1", corr_long = "ar1",
                  algo = "gschur", grid_type = "auto",
                  fpca_fn = NULL, precompute_fpca = "full",
                  fpca_value_col = "resid",
                  fpca_method = c("auto","face","sc"),
                  fpca_sparse_thresh = 0.70,
                  fpca_pve = 0.99,
                  fpca_npc = NULL,
                  fpca_center_sc = FALSE,
                  tol = 1e-8,
                  copy_dt = TRUE,
                  # NEW:
                  ensure_order ="setorderv",# c("none", "setorderv"),
                  check_order = FALSE) {

  fpca_method <- match.arg(fpca_method)
  ensure_order <- match.arg(ensure_order)

  if (second.deriv) stop("This version implements term1 only: set second.deriv = FALSE.")

  dd <- data.table::as.data.table(dd)
  if (copy_dt) dd <- data.table::copy(dd)

  corr_fn   <- tolower(corr_fn)
  corr_long <- tolower(corr_long)

  if (!corr_fn %in% c("ar1", "exchangeable", "independent", "fpca"))
    stop("corr_fn must be one of: 'ar1','exchangeable','independent','fpca'")
  if (!corr_long %in% c("ar1", "exchangeable", "independent"))
    stop("corr_long must be one of: 'ar1','exchangeable','independent'")

  # require upstream columns
  if (!("sqrtv" %in% names(dd))) {
    if (!("v" %in% names(dd))) stop("Need column 'v' or 'sqrtv' in dd.")
    dd[, sqrtv := sqrt(v)]
  }
  if (!("muprime" %in% names(dd))) stop("Need column 'muprime' in dd (compute upstream).")

  # Ensure rho columns exist when needed
  if (corr_long != "independent" && !("rho_long" %in% names(dd))) {
    stop("Need rho_long column in dd (use corr.est upstream).")
  }
  if (corr_fn %in% c("ar1","exchangeable") && !("rho_fn" %in% names(dd))) {
    stop("Need rho_fn column in dd (use corr.est upstream).")
  }

  # ------------------------------------------------------------
  # Optional canonical sort ONCE (no per-group order())
  # ------------------------------------------------------------
  if (ensure_order == "setorderv") {
    data.table::setorderv(dd, c(cname_, index_long, index_fn))
  } else if (isTRUE(check_order)) {
    # Cheap(ish) sanity: within each cluster, (index_long,index_fn) must be nondecreasing lexicographically
    bad <- dd[, {
      il  <- get(index_long)
      ifn <- get(index_fn)
      pil  <- data.table::shift(il)
      pifn <- data.table::shift(ifn)
      any(!is.na(pil) & (il < pil | (il == pil & ifn < pifn)))
    }, by = .(cid = get(cname_))]
    if (any(bad$V1)) {
      stop("dd is not in canonical order within at least one cluster.\n",
           "Expected ordering by (", cname_, ", ", index_long, ", ", index_fn, ").\n",
           "Either sort once upstream, or call .getW(..., ensure_order='setorderv').")
    }
  }

  # ---- FPCA patch: compute if needed ----
  if (corr_fn == "fpca" && is.null(fpca_fn)) {
    fpca_fn <- fgee_fpca_estimate(
      dx = dd,
      cname_ = cname_,
      index_fn = index_fn,
      index_long = index_long,
      value_col = fpca_value_col,
      method = fpca_method,
      sparse_thresh = fpca_sparse_thresh,
      pve = fpca_pve,
      npc = fpca_npc,
      var = TRUE,
      center_sc = fpca_center_sc
    )
  }

  # standardized internals
  dd[, .cid := get(..cname_)]
  dd[, .idx_fn := get(..index_fn)]
  dd[, .idx_long := get(..index_long)]
  dd[, .scl2 := muprime / sqrtv]

  # -------------------------------------------------------------------------
  # Independent x Independent (FAST)
  # -------------------------------------------------------------------------
  if (corr_fn == "independent" && corr_long == "independent") {

    out <- dd[, {
      X <- as.matrix(.SD)
      X <- X * .scl2
      G <- crossprod(X)
      .(W = list(G))
    }, .SDcols = namesd, by = .(.cid)]

    wi <- out$W
    names(wi) <- as.character(out$.cid)

    rm_cols <- intersect(names(dd), c(".cid",".idx_fn",".idx_long",".scl2"))
    dd[, (rm_cols) := NULL]
    return(wi)
  }

  use_kronecker <- (corr_fn != "independent" && corr_long != "independent")

  # -------------------------------------------------------------------------
  # FPCA inverse apply (functional direction)
  # -------------------------------------------------------------------------
  use_precomputed_fpca <- FALSE
  use_grouped_fpca <- FALSE
  apply_fpca_inv_mat <- NULL

  if (corr_fn == "fpca") {

    Phi_full <- as.matrix(fpca_fn$efunctions)
    eigenvalues <- as.numeric(fpca_fn$evalues)
    sigma2_fpca <- max(as.numeric(fpca_fn$sigma2), tol)
    argvals_fpca <- fpca_fn$argvals

    keep <- which(is.finite(eigenvalues) & eigenvalues > 0)
    if (!length(keep)) stop("FPCA eigenvalues are not positive/finite.")
    eigenvalues <- eigenvalues[keep]
    Phi_full <- Phi_full[, keep, drop = FALSE]
    Lambda_inv <- diag(1 / eigenvalues)

    map_idx <- function(idx_vals) {
      pos <- match(idx_vals, argvals_fpca)
      if (anyNA(pos)) stop("FPCA: some index_fn values not in fpca_fn$argvals.")
      pos
    }

    compute_Cinv <- function(grid_idx_vals) {
      pos <- map_idx(grid_idx_vals)
      Phi_grid <- Phi_full[pos, , drop = FALSE]
      M <- Lambda_inv + (1 / sigma2_fpca) * crossprod(Phi_grid)

      E <- eigen(M, symmetric = TRUE)
      tol_e <- max(dim(M)) * .Machine$double.eps * max(E$values)
      inv_vals <- ifelse(E$values > tol_e, 1 / E$values, 0)
      M_inv <- E$vectors %*% (inv_vals * t(E$vectors))

      diag(length(pos)) / sigma2_fpca -
        (Phi_grid %*% M_inv %*% t(Phi_grid)) / (sigma2_fpca^2)
    }

    C_inv_full <- NULL
    C_inv_lookup <- NULL

    if (precompute_fpca == "full") {
      grids <- dd[, .(grid_id = paste(sort(unique(.idx_fn)), collapse = ",")),
                  by = .(.cid, .idx_long)]
      if (data.table::uniqueN(grids$grid_id) == 1) {
        idx_full <- sort(unique(dd$.idx_fn))
        C_inv_full <- compute_Cinv(idx_full)
        use_precomputed_fpca <- TRUE
      }
    } else if (precompute_fpca == "grouped") {
      dd[, grid_sig := paste(sort(.idx_fn), collapse = ","), by = .(.cid, .idx_long)]
      unique_sigs <- unique(dd$grid_sig)
      C_inv_lookup <- lapply(unique_sigs, function(sig) {
        idx <- strsplit(sig, ",")[[1]]
        compute_Cinv(idx)
      })
      names(C_inv_lookup) <- unique_sigs
      use_grouped_fpca <- TRUE
    }

    apply_Cidx_inv_mat <- function(Mrhs, idx_vals) {
      pos <- map_idx(idx_vals)
      Phi_idx <- Phi_full[pos, , drop = FALSE]
      Mmid <- Lambda_inv + (1 / sigma2_fpca) * crossprod(Phi_idx)

      E <- eigen(Mmid, symmetric = TRUE)
      tol_e <- max(dim(Mmid)) * .Machine$double.eps * max(E$values)
      inv_vals <- ifelse(E$values > tol_e, 1 / E$values, 0)
      Mmid_inv <- E$vectors %*% (inv_vals * t(E$vectors))

      Mrhs / sigma2_fpca -
        (Phi_idx %*% (Mmid_inv %*% crossprod(Phi_idx, Mrhs))) / (sigma2_fpca^2)
    }

    apply_fpca_inv_mat <- function(Mrhs, idx_vals, grid_sig = NULL) {
      if (use_precomputed_fpca) return(C_inv_full %*% Mrhs)
      if (use_grouped_fpca) {
        if (is.null(grid_sig) || !nzchar(grid_sig)) grid_sig <- paste(sort(idx_vals), collapse = ",")
        Cinv <- C_inv_lookup[[grid_sig]]
        if (is.null(Cinv)) stop("FPCA grouped: grid_sig not found: ", grid_sig)
        return(Cinv %*% Mrhs)
      }
      apply_Cidx_inv_mat(Mrhs, idx_vals = idx_vals)
    }
  }

  # -------------------------------------------------------------------------
  # 1D correlation inverse apply for MATRIX RHS (n x q)
  # -------------------------------------------------------------------------
  solve_corr_mat <- function(Mrhs, corr, rho, times, grid_type, algo, tol) {

    was_vec <- is.null(dim(Mrhs))
    if (was_vec) Mrhs <- matrix(Mrhs, ncol = 1)

    n <- nrow(Mrhs)
    if (n <= 1) return(if (was_vec) as.vector(Mrhs) else Mrhs)

    rho <- as.numeric(rho)
    if (length(rho) != 1L || !is.finite(rho[1])) stop("rho must be a finite scalar.")
    rho <- rho[1]

    if (corr == "exchangeable") {
      acf <- c(1, rep(rho, n - 1))
      Toep <- SuperGauss::Toeplitz$new(N = n, acf = acf)
      out <- Toep$solve(Mrhs, method = algo, tol = tol)
      return(if (was_vec) as.vector(out) else out)
    }

    is_regular <- TRUE
    if (grid_type == "irregular") is_regular <- FALSE
    if (grid_type == "auto" && n > 1) {
      d <- diff(times)
      is_regular <- (max(d) - min(d)) < 1e-10
    }

    if (is_regular) {
      acf <- rho^(0:(n - 1))
      Toep <- SuperGauss::Toeplitz$new(N = n, acf = acf)
      out <- Toep$solve(Mrhs, method = algo, tol = tol)
      return(if (was_vec) as.vector(out) else out)
    } else {
      if (!requireNamespace("irregulAR1", quietly = TRUE)) {
        stop("Package 'irregulAR1' required for grid_type='irregular'.")
      }
      Q <- (1 / (1 - rho^2)) * irregulAR1::ar1_prec_irregular(sigma = 1, times = times, rho = rho)
      out <- as.matrix(Q %*% Mrhs)
      return(if (was_vec) as.vector(out) else out)
    }
  }

  # -------------------------------------------------------------------------
  # Define grouping variables (no per-group order() anymore)
  # -------------------------------------------------------------------------
  if (!use_kronecker) {
    if (corr_fn == "independent" && corr_long != "independent") {
      dd[, .index_vec := .idx_fn]
      dd[, .ar_idx := .idx_long]
    } else if (corr_fn != "independent" && corr_long == "independent") {
      dd[, .index_vec := .idx_long]
      dd[, .ar_idx := .idx_fn]
    } else {
      stop("Internal: non-kronecker but neither direction is independent. This should not happen.")
    }
  }

  pp <- length(namesd)
  clusts <- unique(dd$.cid)

  # -------------------------------------------------------------------------
  # TERM 1: X^T V^{-1} X  with X scaled by scl2 = muprime/sqrtv
  # ASSUMES dd is already in canonical order inside each group
  # -------------------------------------------------------------------------
  if (use_kronecker) {

    term1 <- dd[, {

      X_raw <- as.matrix(.SD)
      X_mat <- X_raw * .scl2

      n_long <- data.table::uniqueN(.idx_long)
      n_func <- data.table::uniqueN(.idx_fn)

      n_cell <- data.table::uniqueN(data.table::data.table(.idx_long, .idx_fn)) # unique (long,fn)

      if (n_cell != n_long * n_func) {
        stop("Kronecker path requires a complete tensor-product grid within each cluster: ",
             "all combinations of (", index_long, ", ", index_fn, ") must be present. ",
             "In this cluster, unique cells=", n_cell,
             " but n_long*n_func=", n_long * n_func,
             ". This usually means missing functional grid points at some longitudinal indices.")
      }

      if (.N != n_cell) {
        stop("Kronecker path requires exactly one row per (", index_long, ", ", index_fn, ") cell within each cluster. ",
             "Found duplicated cells: total rows .N=", .N, " but unique cells=", n_cell,
             ". This usually means your longitudinal index does not uniquely identify a visit/curve (or you have replicated measurements).")
      }
      times_fn   <- sort(unique(.idx_fn))
      times_long <- sort(unique(.idx_long))

      if (corr_fn == "fpca") {
        Ainv <- function(M) apply_fpca_inv_mat(M, idx_vals = times_fn, grid_sig = NULL)
      } else {
        rho_fn_val <- unique(rho_fn)
        if (length(rho_fn_val) != 1) stop("rho_fn must be constant within cluster for Kronecker path.")
        Ainv <- function(M) solve_corr_mat(M, corr = corr_fn, rho = rho_fn_val,
                                           times = times_fn, grid_type = grid_type, algo = algo, tol = tol)
      }

      rho_long_val <- unique(rho_long)
      if (length(rho_long_val) != 1) stop("rho_long must be constant within cluster for Kronecker path.")
      Binv <- function(M) solve_corr_mat(M, corr = corr_long, rho = rho_long_val,
                                         times = times_long, grid_type = grid_type, algo = algo, tol = tol)

      Vinv_cols <- lapply(seq_len(ncol(X_mat)), function(k) {
        Xk <- matrix(X_mat[, k], nrow = n_func, ncol = n_long, byrow = FALSE)
        AX <- Ainv(Xk)
        Z  <- Binv(t(AX))
        as.vector(t(Z))
      })

      Vinv_X <- do.call(cbind, Vinv_cols)
      G <- crossprod(X_mat, Vinv_X)
      vec <- as.vector(G)
      .(m_idx = seq_along(vec), T1 = vec)
    }, .SDcols = namesd, by = .(.cid)]

    term1_sum <- term1

  } else {

    if (corr_fn == "fpca" && corr_long == "independent") {

      term1 <- dd[, {

        X_raw <- as.matrix(.SD)
        X_mat <- X_raw * .scl2
        idx_vals <- .ar_idx

        grid_sig_val <- if ("grid_sig" %in% names(dd)) unique(grid_sig) else NULL
        Vinv_X <- apply_fpca_inv_mat(X_mat, idx_vals = idx_vals, grid_sig = grid_sig_val)

        G <- crossprod(X_mat, Vinv_X)
        vec <- as.vector(G)
        .(m_idx = seq_along(vec), T1 = vec)
      }, .SDcols = namesd, by = .(.cid, .index_vec)]

      term1_sum <- term1[, .(T1 = sum(T1)), by = .(.cid, m_idx)]

    } else {

      if (corr_fn == "independent") {
        rho_var <- "rho_long"
        corr_type <- corr_long
      } else {
        rho_var <- "rho_fn"
        corr_type <- corr_fn
      }

      term1 <- dd[, {

        X_raw <- as.matrix(.SD)
        X_mat <- X_raw * .scl2
        times <- .ar_idx

        rho_val <- unique(get(rho_var))
        if (length(rho_val) != 1) stop(rho_var, " must be constant within (cluster, index_vec).")

        Vinv_X <- solve_corr_mat(X_mat, corr = corr_type, rho = rho_val,
                                 times = times, grid_type = grid_type, algo = algo, tol = tol)

        G <- crossprod(X_mat, Vinv_X)
        vec <- as.vector(G)
        .(m_idx = seq_along(vec), T1 = vec)
      }, .SDcols = namesd, by = .(.cid, .index_vec)]

      term1_sum <- term1[, .(T1 = sum(T1)), by = .(.cid, m_idx)]
    }
  }

  out <- lapply(clusts, function(cc) {
    vec <- term1_sum[.cid == cc][order(m_idx), T1]
    matrix(vec, byrow = TRUE, nrow = pp, ncol = pp) # took out negative
  })
  names(out) <- as.character(clusts)

  rm_cols <- intersect(names(dd), c(".cid",".idx_fn",".idx_long",".scl2",".index_vec",".ar_idx","grid_sig"))
  dd[, (rm_cols) := NULL]

  out
}




# ============================================================================
# .getD()
# ============================================================================
#' @keywords internal
#' @noRd
.getD <- function(dd, namesd, cname_,
                  index_fn = "yindex.vec", index_long = "time",
                  corr_fn = "ar1", corr_long = "ar1",
                  algo = "gschur", grid_type = "auto",
                  fpca_fn = NULL, precompute_fpca = "full",
                  fpca_value_col = "resid",
                  fpca_method = c("auto","face","sc"),
                  fpca_sparse_thresh = 0.70,
                  fpca_pve = 0.99,
                  fpca_npc = NULL,
                  fpca_center_sc = FALSE,
                  tol = 1e-8,
                  copy_dt = TRUE,
                  resid_col = "resid",
                  # NEW:
                  ensure_order = "setorderv",#c("none", "setorderv"),
                  check_order = FALSE) {

  fpca_method <- match.arg(fpca_method)
  ensure_order <- match.arg(ensure_order)

  dd <- data.table::as.data.table(dd)
  if (copy_dt) dd <- data.table::copy(dd)

  corr_fn   <- tolower(corr_fn)
  corr_long <- tolower(corr_long)

  if (!corr_fn %in% c("ar1", "exchangeable", "independent", "fpca"))
    stop("corr_fn must be one of: 'ar1','exchangeable','independent','fpca'")
  if (!corr_long %in% c("ar1", "exchangeable", "independent"))
    stop("corr_long must be one of: 'ar1','exchangeable','independent'")

  if (!(resid_col %in% names(dd))) stop("Need resid_col='", resid_col, "' in dd (compute upstream).")
  if (!("sqrtv" %in% names(dd))) {
    if (!("v" %in% names(dd))) stop("Need column 'v' or 'sqrtv' in dd.")
    dd[, sqrtv := sqrt(v)]
  }
  if (!("muprime" %in% names(dd))) stop("Need column 'muprime' in dd (compute upstream).")

  # Ensure rho columns exist when needed
  if (corr_long != "independent" && !("rho_long" %in% names(dd))) {
    stop("Need rho_long column in dd (use corr.est upstream).")
  }
  if (corr_fn %in% c("ar1","exchangeable") && !("rho_fn" %in% names(dd))) {
    stop("Need rho_fn column in dd (use corr.est upstream).")
  }

  # ------------------------------------------------------------
  # Optional canonical sort ONCE (no per-group order())
  # ------------------------------------------------------------
  if (ensure_order == "setorderv") {
    data.table::setorderv(dd, c(cname_, index_long, index_fn))
  } else if (isTRUE(check_order)) {
    bad <- dd[, {
      il  <- get(index_long)
      ifn <- get(index_fn)
      pil  <- data.table::shift(il)
      pifn <- data.table::shift(ifn)
      any(!is.na(pil) & (il < pil | (il == pil & ifn < pifn)))
    }, by = .(cid = get(cname_))]
    if (any(bad$V1)) {
      stop("dd is not in canonical order within at least one cluster.\n",
           "Expected ordering by (", cname_, ", ", index_long, ", ", index_fn, ").\n",
           "Either sort once upstream, or call .getD(..., ensure_order='setorderv').")
    }
  }

  if (corr_fn == "fpca" && is.null(fpca_fn)) {
    fpca_fn <- fgee_fpca_estimate(
      dx = dd,
      cname_ = cname_,
      index_fn = index_fn,
      index_long = index_long,
      value_col = fpca_value_col,
      method = fpca_method,
      sparse_thresh = fpca_sparse_thresh,
      pve = fpca_pve,
      npc = fpca_npc,
      var = TRUE,
      center_sc = fpca_center_sc
    )
  }

  # standardized internals
  dd[, .cid := get(..cname_)]
  dd[, .idx_fn := get(..index_fn)]
  dd[, .idx_long := get(..index_long)]
  dd[, .resid := get(resid_col)]

  # -------------------------------------------------------------------------
  # Independent x Independent (FAST)
  # -------------------------------------------------------------------------
  if (corr_fn == "independent" && corr_long == "independent") {

    dd[, .vinv_r := .resid / sqrtv]  # (Y-mu)/v

    out <- dd[, {
      X <- as.matrix(.SD)
      X <- X * muprime
      vec <- as.numeric(crossprod(X, .vinv_r))
      .(D = list(vec))
    }, .SDcols = namesd, by = .(.cid)]

    di <- out$D
    names(di) <- as.character(out$.cid)

    rm_cols <- intersect(names(dd), c(".cid",".idx_fn",".idx_long",".resid",".vinv_r"))
    dd[, (rm_cols) := NULL]
    return(di)
  }

  use_kronecker <- (corr_fn != "independent" && corr_long != "independent")

  # -------------------------------------------------------------------------
  # FPCA inverse apply (functional direction)
  # -------------------------------------------------------------------------
  use_precomputed_fpca <- FALSE
  use_grouped_fpca <- FALSE
  apply_fpca_inv_mat <- NULL

  if (corr_fn == "fpca") {

    Phi_full <- as.matrix(fpca_fn$efunctions)
    eigenvalues <- as.numeric(fpca_fn$evalues)
    sigma2_fpca <- max(as.numeric(fpca_fn$sigma2), tol)
    argvals_fpca <- fpca_fn$argvals

    keep <- which(is.finite(eigenvalues) & eigenvalues > 0)
    if (!length(keep)) stop("FPCA eigenvalues are not positive/finite.")
    eigenvalues <- eigenvalues[keep]
    Phi_full <- Phi_full[, keep, drop = FALSE]
    Lambda_inv <- diag(1 / eigenvalues)

    map_idx <- function(idx_vals) {
      pos <- match(idx_vals, argvals_fpca)
      if (anyNA(pos)) stop("FPCA: some index_fn values not in fpca_fn$argvals.")
      pos
    }

    compute_Cinv <- function(grid_idx_vals) {
      pos <- map_idx(grid_idx_vals)
      Phi_grid <- Phi_full[pos, , drop = FALSE]
      M <- Lambda_inv + (1 / sigma2_fpca) * crossprod(Phi_grid)

      E <- eigen(M, symmetric = TRUE)
      tol_e <- max(dim(M)) * .Machine$double.eps * max(E$values)
      inv_vals <- ifelse(E$values > tol_e, 1 / E$values, 0)
      M_inv <- E$vectors %*% (inv_vals * t(E$vectors))

      diag(length(pos)) / sigma2_fpca -
        (Phi_grid %*% M_inv %*% t(Phi_grid)) / (sigma2_fpca^2)
    }

    C_inv_full <- NULL
    C_inv_lookup <- NULL

    if (precompute_fpca == "full") {
      grids <- dd[, .(grid_id = paste(sort(unique(.idx_fn)), collapse = ",")),
                  by = .(.cid, .idx_long)]
      if (data.table::uniqueN(grids$grid_id) == 1) {
        idx_full <- sort(unique(dd$.idx_fn))
        C_inv_full <- compute_Cinv(idx_full)
        use_precomputed_fpca <- TRUE
      }
    } else if (precompute_fpca == "grouped") {
      dd[, grid_sig := paste(sort(.idx_fn), collapse = ","), by = .(.cid, .idx_long)]
      unique_sigs <- unique(dd$grid_sig)
      C_inv_lookup <- lapply(unique_sigs, function(sig) {
        idx <- strsplit(sig, ",")[[1]]
        compute_Cinv(idx)
      })
      names(C_inv_lookup) <- unique_sigs
      use_grouped_fpca <- TRUE
    }

    apply_Cidx_inv_mat <- function(Mrhs, idx_vals) {
      pos <- map_idx(idx_vals)
      Phi_idx <- Phi_full[pos, , drop = FALSE]
      Mmid <- Lambda_inv + (1 / sigma2_fpca) * crossprod(Phi_idx)

      E <- eigen(Mmid, symmetric = TRUE)
      tol_e <- max(dim(Mmid)) * .Machine$double.eps * max(E$values)
      inv_vals <- ifelse(E$values > tol_e, 1 / E$values, 0)
      Mmid_inv <- E$vectors %*% (inv_vals * t(E$vectors))

      Mrhs / sigma2_fpca -
        (Phi_idx %*% (Mmid_inv %*% crossprod(Phi_idx, Mrhs))) / (sigma2_fpca^2)
    }

    apply_fpca_inv_mat <- function(Mrhs, idx_vals, grid_sig = NULL) {
      if (use_precomputed_fpca) return(C_inv_full %*% Mrhs)
      if (use_grouped_fpca) {
        if (is.null(grid_sig) || !nzchar(grid_sig)) grid_sig <- paste(sort(idx_vals), collapse = ",")
        Cinv <- C_inv_lookup[[grid_sig]]
        if (is.null(Cinv)) stop("FPCA grouped: grid_sig not found: ", grid_sig)
        return(Cinv %*% Mrhs)
      }
      apply_Cidx_inv_mat(Mrhs, idx_vals = idx_vals)
    }
  }

  # -------------------------------------------------------------------------
  # 1D correlation inverse apply
  # -------------------------------------------------------------------------
  solve_corr_mat <- function(Mrhs, corr, rho, times, grid_type, algo, tol) {
    was_vec <- is.null(dim(Mrhs))
    if (was_vec) Mrhs <- matrix(Mrhs, ncol = 1)

    n <- nrow(Mrhs)
    if (n <= 1) return(if (was_vec) as.vector(Mrhs) else Mrhs)

    rho <- as.numeric(rho)
    if (length(rho) != 1L || !is.finite(rho[1])) stop("rho must be a finite scalar.")
    rho <- rho[1]

    if (corr == "exchangeable") {
      acf <- c(1, rep(rho, n - 1))
      Toep <- SuperGauss::Toeplitz$new(N = n, acf = acf)
      out <- Toep$solve(Mrhs, method = algo, tol = tol)
      return(if (was_vec) as.vector(out) else out)
    }

    is_regular <- TRUE
    if (grid_type == "irregular") is_regular <- FALSE
    if (grid_type == "auto" && n > 1) {
      d <- diff(times)
      is_regular <- (max(d) - min(d)) < 1e-10
    }

    if (is_regular) {
      acf <- rho^(0:(n - 1))
      Toep <- SuperGauss::Toeplitz$new(N = n, acf = acf)
      out <- Toep$solve(Mrhs, method = algo, tol = tol)
      return(if (was_vec) as.vector(out) else out)
    } else {
      if (!requireNamespace("irregulAR1", quietly = TRUE)) {
        stop("Package 'irregulAR1' required for grid_type='irregular'.")
      }
      Q <- (1 / (1 - rho^2)) * irregulAR1::ar1_prec_irregular(sigma = 1, times = times, rho = rho)
      out <- as.matrix(Q %*% Mrhs)
      return(if (was_vec) as.vector(out) else out)
    }
  }

  # -------------------------------------------------------------------------
  # Define grouping variables (no per-group order() anymore)
  # -------------------------------------------------------------------------
  if (!use_kronecker) {
    if (corr_fn == "independent" && corr_long != "independent") {
      dd[, .index_vec := .idx_fn]
      dd[, .ar_idx := .idx_long]
    } else if (corr_fn != "independent" && corr_long == "independent") {
      dd[, .index_vec := .idx_long]
      dd[, .ar_idx := .idx_fn]
    } else {
      stop("Internal: non-kronecker but neither direction is independent. This should not happen.")
    }
  }

  clusts <- unique(dd$.cid)

  # -------------------------------------------------------------------------
  # CORE: D^T V^{-1} r   with D = diag(muprime) X
  # resid already includes 1/sqrtv; apply R^{-1}; then /sqrtv
  # ASSUMES dd is already in canonical order inside each group
  # -------------------------------------------------------------------------
  if (use_kronecker) {

    res <- dd[, {

      r0 <- .resid
      sv <- sqrtv

      n_long <- data.table::uniqueN(.idx_long)
      n_func <- data.table::uniqueN(.idx_fn)

      n_cell <- data.table::uniqueN(data.table::data.table(.idx_long, .idx_fn)) # unique (long,fn)

      if (n_cell != n_long * n_func) {
        stop("Kronecker path requires a complete tensor-product grid within each cluster: ",
             "all combinations of (", index_long, ", ", index_fn, ") must be present. ",
             "In this cluster, unique cells=", n_cell,
             " but n_long*n_func=", n_long * n_func,
             ". This usually means missing functional grid points at some longitudinal indices.")
      }

      if (.N != n_cell) {
        stop("Kronecker path requires exactly one row per (", index_long, ", ", index_fn, ") cell within each cluster. ",
             "Found duplicated cells: total rows .N=", .N, " but unique cells=", n_cell,
             ". This usually means your longitudinal index does not uniquely identify a visit/curve (or you have replicated measurements).")
      }

      times_fn   <- sort(unique(.idx_fn))
      times_long <- sort(unique(.idx_long))

      if (corr_fn == "fpca") {
        Ainv <- function(M) apply_fpca_inv_mat(M, idx_vals = times_fn, grid_sig = NULL)
      } else {
        rho_fn_val <- unique(rho_fn)
        if (length(rho_fn_val) != 1) stop("rho_fn must be constant within cluster in Kronecker path.")
        Ainv <- function(M) solve_corr_mat(M, corr = corr_fn, rho = rho_fn_val,
                                           times = times_fn, grid_type = grid_type, algo = algo, tol = tol)
      }

      rho_long_val <- unique(rho_long)
      if (length(rho_long_val) != 1) stop("rho_long must be constant within cluster in Kronecker path.")
      Binv <- function(M) solve_corr_mat(M, corr = corr_long, rho = rho_long_val,
                                         times = times_long, grid_type = grid_type, algo = algo, tol = tol)

      Rhs <- matrix(r0, nrow = n_func, ncol = n_long, byrow = FALSE)
      A_Rhs <- Ainv(Rhs)
      Z <- Binv(t(A_Rhs))
      rRinv <- as.vector(t(Z))

      Vinv_r <- rRinv / sv

      X_raw <- as.matrix(.SD)
      X_mat <- X_raw * muprime
      vec <- as.numeric(crossprod(X_mat, Vinv_r))
      .(m_idx = seq_along(vec), V1 = vec)
    }, .SDcols = namesd, by = .(.cid)]

    out <- lapply(clusts, function(ii) as.numeric(res[.cid == ii][order(m_idx), V1]))
    names(out) <- as.character(clusts)

    rm_cols <- intersect(names(dd), c(".cid",".idx_fn",".idx_long",".resid",".index_vec",".ar_idx","grid_sig"))
    dd[, (rm_cols) := NULL]
    return(out)

  } else {

    if (corr_fn == "fpca" && corr_long == "independent") {

      res <- dd[, {

        r0 <- .resid
        sv <- sqrtv
        idx_vals <- .ar_idx

        rRinv <- as.vector(apply_fpca_inv_mat(matrix(r0, ncol = 1),
                                              idx_vals = idx_vals,
                                              grid_sig = if ("grid_sig" %in% names(dd)) unique(grid_sig) else NULL))
        Vinv_r <- rRinv / sv

        X_raw <- as.matrix(.SD)
        X_mat <- X_raw * muprime
        vec <- as.numeric(crossprod(X_mat, Vinv_r))
        .(m_idx = seq_along(vec), V1 = vec)
      }, .SDcols = namesd, by = .(.cid, .index_vec)][
        , .(V1 = sum(V1)), by = .(.cid, m_idx)
      ]

      out <- lapply(clusts, function(ii) as.numeric(res[.cid == ii][order(m_idx), V1]))
      names(out) <- as.character(clusts)

      rm_cols <- intersect(names(dd), c(".cid",".idx_fn",".idx_long",".resid",".index_vec",".ar_idx","grid_sig"))
      dd[, (rm_cols) := NULL]
      return(out)

    } else {

      if (corr_fn == "independent") {
        rho_var <- "rho_long"
        corr_type <- corr_long
      } else {
        rho_var <- "rho_fn"
        corr_type <- corr_fn
      }

      res <- dd[, {

        r0 <- .resid
        sv <- sqrtv
        times <- .ar_idx

        rho_val <- unique(get(rho_var))
        if (length(rho_val) != 1L || !is.finite(rho_val)) {
          stop(sprintf("%s must be constant and finite within each (cluster, index_vec) block.", rho_var))
        }

        rRinv <- solve_corr_mat(r0, corr = corr_type, rho = rho_val,
                                times = times, grid_type = grid_type, algo = algo, tol = tol)

        Vinv_r <- rRinv / sv

        X_raw <- as.matrix(.SD)
        X_mat <- X_raw * muprime
        vec <- as.numeric(crossprod(X_mat, Vinv_r))
        .(m_idx = seq_along(vec), V1 = vec)
      }, .SDcols = namesd, by = .(.cid, .index_vec)][
        , .(V1 = sum(V1)), by = .(.cid, m_idx)
      ]

      out <- lapply(clusts, function(ii) as.numeric(res[.cid == ii][order(m_idx), V1]))
      names(out) <- as.character(clusts)

      rm_cols <- intersect(names(dd), c(".cid",".idx_fn",".idx_long",".resid",".index_vec",".ar_idx","grid_sig"))
      dd[, (rm_cols) := NULL]
      return(out)
    }
  }
}



# Internal functions
#' @keywords internal
#' @noRd
W.estimate <- function(dx,
                       namesd,
                       cname_,
                       corr_fn = "independent",
                       corr_long = "independent",
                       index_fn = "yindex.vec",
                       index_long = "time",
                       algo = "gschur",
                       grid_type = "auto",
                       fpca_fn = NULL,
                       precompute_fpca = "full",
                       tol = 1e-8,
                       copy_dt = TRUE,
                       fallback_algo = c("pcg"),
                       id.vec = NULL) {

  dd <- data.table::as.data.table(dx)
  if (copy_dt) dd <- data.table::copy(dd)

  wi <- .getW(dd, namesd = namesd, cname_ = cname_,
              index_fn = index_fn, index_long = index_long,
              corr_fn = corr_fn, corr_long = corr_long,
              algo = algo, grid_type = grid_type,
              fpca_fn = fpca_fn, precompute_fpca = precompute_fpca,
              tol = tol, copy_dt = FALSE)

  bad <- which(vapply(wi, function(M) any(!is.finite(M)), logical(1)))
  if (length(bad) > 0 && length(fallback_algo) > 0) {

    for (a2 in fallback_algo) {
      wi2 <- .getW(dd, namesd = namesd, cname_ = cname_,
                   index_fn = index_fn, index_long = index_long,
                   corr_fn = corr_fn, corr_long = corr_long,
                   algo = a2, grid_type = grid_type,
                   fpca_fn = fpca_fn, precompute_fpca = precompute_fpca,
                   tol = tol, copy_dt = FALSE)

      bad2 <- which(vapply(wi2, function(M) any(!is.finite(M)), logical(1)))
      if (length(bad2) == 0) {
        wi <- wi2
        algo <- a2
        break
      }
    }
  }

  bad <- which(vapply(wi, function(M) any(!is.finite(M)), logical(1)))
  if (length(bad) > 0) {
    if (!is.null(id.vec) && length(id.vec) >= max(bad)) {
      print(paste("These clusters/IDs yield errors (please remove):", paste(id.vec[bad], collapse = ",")))
    }
    stop("W.estimate returned non-finite values for some clusters.")
  }

  attr(wi, "algo") <- algo
  wi
}


#' @keywords internal
#' @noRd
D.estimate <- function(dx,
                       namesd,
                       cname_,
                       corr_fn = "independent",
                       corr_long = "independent",
                       index_fn = "yindex.vec",
                       index_long = "time",
                       algo = "gschur",
                       grid_type = "auto",
                       fpca_fn = NULL,
                       precompute_fpca = "full",
                       tol = 1e-8,
                       copy_dt = TRUE,
                       fallback_algo = c("pcg"),
                       id.vec = NULL,
                       resid_col = "resid") {

  dd <- data.table::as.data.table(dx)
  if (copy_dt) dd <- data.table::copy(dd)

  di <- .getD(dd, namesd = namesd, cname_ = cname_,
              index_fn = index_fn, index_long = index_long,
              corr_fn = corr_fn, corr_long = corr_long,
              algo = algo, grid_type = grid_type,
              fpca_fn = fpca_fn, precompute_fpca = precompute_fpca,
              tol = tol, copy_dt = FALSE, resid_col = resid_col)

  bad <- which(vapply(di, function(v) any(!is.finite(v)), logical(1)))
  if (length(bad) > 0 && length(fallback_algo) > 0) {

    for (a2 in fallback_algo) {
      di2 <- .getD(dd, namesd = namesd, cname_ = cname_,
                   index_fn = index_fn, index_long = index_long,
                   corr_fn = corr_fn, corr_long = corr_long,
                   algo = a2, grid_type = grid_type,
                   fpca_fn = fpca_fn, precompute_fpca = precompute_fpca,
                   tol = tol, copy_dt = FALSE, resid_col = resid_col)

      bad2 <- which(vapply(di2, function(v) any(!is.finite(v)), logical(1)))
      if (length(bad2) == 0) {
        di <- di2
        algo <- a2
        break
      }
    }
  }

  bad <- which(vapply(di, function(v) any(!is.finite(v)), logical(1)))
  if (length(bad) > 0) {
    if (!is.null(id.vec) && length(id.vec) >= max(bad)) {
      print(paste("These clusters/IDs yield errors (please remove):", paste(id.vec[bad], collapse = ",")))
    }
    stop("D.estimate returned non-finite values for some clusters.")
  }

  attr(di, "algo") <- algo
  di
}
