#' Estimate working correlations in longitudinal and/or functional directions.
#' Produces rho_long and rho_fn columns needed by .getW()/.getD().
#'
#' Rules:
#' - If corr_long != independent and corr_fn == independent: rho_long can vary over index_fn.
#' - If corr_fn != independent and corr_long == independent: rho_fn can vary over index_long.
#' - If BOTH directions non-independent (Kronecker/separable): rho_long and rho_fn are pooled
#'   to cluster-level scalars (constant across the other dimension).
#' - If corr_fn == "fpca": rho_fn is not used (set NA); fpca computed optionally.
#'
#' Estimate working correlations in longitudinal and/or functional directions.
#' Produces rho_long and rho_fn columns needed by .getW()/.getD().
#'
#' IMPORTANT: This version POOLS across clusters ALWAYS, so rho_* are identical
#' across clusters (either a scalar or a function of index_fn/index_long).
#'
#' @keywords internal
#' @noRd
corr.est <- function(dx,
                     cname_ = "cname_",
                     index_fn = "yindex.vec",
                     index_long = "time",
                     corr_fn = "independent",
                     corr_long = "independent",
                     resid_col = "resid",
                     rho.smooth = FALSE,
                     ar = c("mom", "yw"),
                     glmfit = NULL,
                     clamp = 0.999,
                     copy_dt = TRUE,
                     # FPCA (optional):
                     fpca_fn = NULL,
                     fpca_method = c("auto","face","sc"),
                     fpca_sparse_thresh = 0.70,
                     fpca_pve = 0.99,
                     fpca_npc = NULL,
                     fpca_center_sc = FALSE,
                     verbose = FALSE,
                     # NEW: avoid expensive order() when times already sorted
                     sort_times = c("auto","always","never")) {

  ar <- match.arg(ar)
  fpca_method <- match.arg(fpca_method)
  sort_times <- match.arg(sort_times)

  dt <- data.table::as.data.table(dx)
  if (copy_dt) dt <- data.table::copy(dt)

  corr_fn   <- tolower(corr_fn)
  corr_long <- tolower(corr_long)

  if (!corr_fn %in% c("ar1","exchangeable","independent","fpca"))
    stop("corr_fn must be one of: ar1, exchangeable, independent, fpca")
  if (!corr_long %in% c("ar1","exchangeable","independent"))
    stop("corr_long must be one of: ar1, exchangeable, independent")

  if (!(cname_ %in% names(dt))) stop("Missing cname_ column: ", cname_)
  if (!(index_fn %in% names(dt))) stop("Missing index_fn column: ", index_fn)
  if (!(index_long %in% names(dt))) stop("Missing index_long column: ", index_long)
  if (!(resid_col %in% names(dt))) stop("Missing resid_col column: ", resid_col)

  # internal standard names
  dt[, .cid      := get(..cname_)]
  dt[, .idx_fn   := get(..index_fn)]
  dt[, .idx_long := get(..index_long)]
  dt[, .resid    := get(..resid_col)]

  if (!("rho_long" %in% names(dt))) dt[, rho_long := NA_real_]
  if (!("rho_fn"   %in% names(dt))) dt[, rho_fn   := NA_real_]

  kronecker <- (corr_fn != "independent" && corr_long != "independent")

  # ----- helpers -----
  .clamp_rho <- function(rho, type) {
    # Handle edge cases
    if (length(rho) == 0) return(numeric(0))

    # Replace non-finite values with 0
    rho <- ifelse(!is.finite(rho), 0, rho)

    # Clamp based on type
    if (type == "ar1") {
      pmax(0, pmin(rho, clamp))
    } else {  # exchangeable
      ifelse(abs(rho) > clamp, clamp * sign(rho), rho)
    }
  }

  .ar1_numden <- function(r, t, ar) {
    n <- length(r)
    if (n < 2) return(list(num = NA_real_, den = NA_real_))

    # NEW: only sort if needed (or if forced)
    if (sort_times == "always") {
      ord <- order(t)
      r <- r[ord]
    } else if (sort_times == "auto") {
      if (is.unsorted(t)) {
        ord <- order(t)
        r <- r[ord]
      }
    } else {
      # sort_times == "never": assume caller already provided sorted times
      # (fastest, but can be wrong if unsorted)
    }

    if (ar == "yw") {
      rho_i <- tryCatch({
        out <- Rfast::ar1(r, method = "yw")
        if (is.numeric(out) && length(out) >= 2) as.numeric(out[2]) else NA_real_
      }, error = function(e) NA_real_)
      if (!is.finite(rho_i)) return(list(num = NA_real_, den = NA_real_))
      return(list(num = rho_i, den = 1))
    }

    num <- sum(r[-1] * r[-n])
    den <- sum(r[-n]^2)
    if (!is.finite(num) || !is.finite(den) || den <= 0) return(list(num = NA_real_, den = NA_real_))
    list(num = num, den = den)
  }

  .ex_numden <- function(r) {
    n <- length(r)
    if (n < 2) return(list(num = NA_real_, den = NA_real_))
    s1 <- sum(r)
    s2 <- sum(r^2)
    num <- (s1^2 - s2)   # = 2*sum_{i<j} r_i r_j
    den <- n * (n - 1)
    if (!is.finite(num) || !is.finite(den) || den <= 0) return(list(num = NA_real_, den = NA_real_))
    list(num = num, den = den)
  }

  # ------------------------------------------------------------
  # rho_long
  # ------------------------------------------------------------
  rho_long_tbl <- NULL

  if (corr_long == "independent") {

    dt[, rho_long := 0]

  } else if (kronecker) {

    if (verbose) message("corr.est: Kronecker pooled rho_long")

    if (corr_long == "ar1") {
      st <- dt[, {
        tmp <- .ar1_numden(.resid, .idx_long, ar = ar)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_fn)]
    } else {
      st <- dt[, {
        tmp <- .ex_numden(.resid)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_fn)]
    }

    st <- st[is.finite(num) & is.finite(den) & den > 0]
    rhoL <- if (nrow(st) > 0) sum(st$num) / sum(st$den) else 0
    rhoL <- .clamp_rho(rhoL, corr_long)

    dt[, rho_long := rhoL]
    rho_long_tbl <- data.table::data.table(rho_long = rhoL)

  } else {

    if (verbose) message("corr.est: rho_long(index_fn), pooled across clusters")

    if (corr_long == "ar1") {
      st <- dt[, {
        tmp <- .ar1_numden(.resid, .idx_long, ar = ar)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_fn)]
    } else {
      st <- dt[, {
        tmp <- .ex_numden(.resid)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_fn)]
    }

    st <- st[is.finite(num) & is.finite(den) & den > 0]
    rho_long_tbl <- st[, .(rho_long = sum(num) / sum(den)), by = .(.idx_fn)]
    rho_long_tbl[, rho_long := .clamp_rho(rho_long, corr_long)]

    if (isTRUE(rho.smooth)) {
      if (is.null(glmfit)) stop("rho.smooth=TRUE requires glmfit.")
      if (is.numeric(rho_long_tbl$.idx_fn) && nrow(rho_long_tbl) >= 6) {
        k0 <- min(glmfit$smooth[[1]]$bs.dim + 1,
                  data.table::uniqueN(rho_long_tbl$.idx_fn) - 1)
        if (is.finite(k0) && k0 >= 4) {
          fit <- mgcv::gam(rho_long ~ s(.idx_fn, k = k0), data = rho_long_tbl)
          rho_long_tbl[, rho_long := as.numeric(fit$fitted.values)]
          rho_long_tbl[, rho_long := .clamp_rho(rho_long, corr_long)]
        }
      }
    }

    dt[rho_long_tbl, rho_long := i.rho_long, on = .(.idx_fn)]
  }

  # ------------------------------------------------------------
  # rho_fn  (unless fpca)
  # ------------------------------------------------------------
  rho_fn_tbl <- NULL

  if (corr_fn == "independent") {

    dt[, rho_fn := 0]

  } else if (corr_fn == "fpca") {

    dt[, rho_fn := NA_real_]

  } else if (kronecker) {

    if (verbose) message("corr.est: Kronecker pooled rho_fn")

    if (corr_fn == "ar1") {
      st <- dt[, {
        tmp <- .ar1_numden(.resid, .idx_fn, ar = ar)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_long)]
    } else {
      st <- dt[, {
        tmp <- .ex_numden(.resid)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_long)]
    }

    st <- st[is.finite(num) & is.finite(den) & den > 0]
    rhoF <- if (nrow(st) > 0) sum(st$num) / sum(st$den) else 0
    rhoF <- .clamp_rho(rhoF, corr_fn)

    dt[, rho_fn := rhoF]
    rho_fn_tbl <- data.table::data.table(rho_fn = rhoF)

  } else {

    if (verbose) message("corr.est: rho_fn(index_long), pooled across clusters")

    if (corr_fn == "ar1") {
      st <- dt[, {
        tmp <- .ar1_numden(.resid, .idx_fn, ar = ar)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_long)]
    } else {
      st <- dt[, {
        tmp <- .ex_numden(.resid)
        .(num = tmp$num, den = tmp$den)
      }, by = .(.cid, .idx_long)]
    }

    st <- st[is.finite(num) & is.finite(den) & den > 0]
    rho_fn_tbl <- st[, .(rho_fn = sum(num) / sum(den)), by = .(.idx_long)]
    rho_fn_tbl[, rho_fn := .clamp_rho(rho_fn, corr_fn)]

    if (isTRUE(rho.smooth)) {
      if (is.null(glmfit)) stop("rho.smooth=TRUE requires glmfit.")
      if (is.numeric(rho_fn_tbl$.idx_long) && nrow(rho_fn_tbl) >= 6) {
        k0 <- min(glmfit$smooth[[1]]$bs.dim + 1,
                  data.table::uniqueN(rho_fn_tbl$.idx_long) - 1)
        if (is.finite(k0) && k0 >= 4) {
          fit <- mgcv::gam(rho_fn ~ s(.idx_long, k = k0), data = rho_fn_tbl)
          rho_fn_tbl[, rho_fn := as.numeric(fit$fitted.values)]
          rho_fn_tbl[, rho_fn := .clamp_rho(rho_fn, corr_fn)]
        }
      }
    }

    dt[rho_fn_tbl, rho_fn := i.rho_fn, on = .(.idx_long)]
  }

  # ------------------------------------------------------------
  # FPCA (optional) if corr_fn == fpca
  # ------------------------------------------------------------
  if (corr_fn == "fpca" && is.null(fpca_fn)) {
    if (verbose) message("corr.est: computing FPCA (corr_fn='fpca')")

    # NEW: protect dt from any by-reference modifications inside FPCA helper
    fpca_fn <- fgee_fpca_estimate(
      dx = dt,
      cname_ = ".cid",
      index_fn = ".idx_fn",
      index_long = ".idx_long",
      value_col = ".resid",
      method = fpca_method,
      sparse_thresh = fpca_sparse_thresh,
      pve = fpca_pve,
      npc = fpca_npc,
      var = TRUE,
      center_sc = fpca_center_sc
    )
  }

  dt[, c(".cid", ".idx_fn", ".idx_long", ".resid") := NULL]

  list(
    dx = dt,
    fpca = fpca_fn,
    rho = list(rho_long = rho_long_tbl, rho_fn = rho_fn_tbl),
    info = list(kronecker = kronecker, corr_fn = corr_fn, corr_long = corr_long,
                sort_times = sort_times)
  )
}


