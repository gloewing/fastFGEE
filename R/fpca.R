#' Estimate FPCA object for corr_fn = "fpca" from long-format dx
#'
#' Returns a list with efunctions/evalues/sigma2 suitable for fpca_fn
#' argument in .getW_AR_dt() / .getD_AR_dt().
#'
#' @param dx data.frame/data.table in LONG format with columns:
#'   cluster id (cname_), functional index (index_fn), optional longitudinal index (index_long),
#'   and the residual column (value_col) you want FPCA on (usually standardized residuals).
#' @param cname_ character, cluster column name
#' @param index_fn character, functional-domain index column (e.g. "yindex.vec")
#' @param index_long character or NULL, longitudinal index column (e.g. "time").
#'   If present, curves are defined by (cluster, time). If NULL or missing, curves by cluster only.
#' @param value_col character, column containing the values to FPCA (recommend: scaled residuals "resid")
#' @param method "auto" (default), "face", or "sc"
#' @param sparse_thresh heuristic threshold (median observed fraction) below which use fpca.sc
#' @param center passed to refund FPCA
#' @param pve, npc passed to refund FPCA (see docs)
#' @param face_args list of extra args for fpca.face
#' @param sc_args list of extra args for fpca.sc
#' @param keep_fit keep full fpca.* fit object
#' @param verbose message progress / decision
#'
#' @return list with efunctions, evalues, sigma2, argvals, mu, npc, pve, method_used
# -----------------------------------------------------------------------------
# FPCA estimator for functional-direction correlation
# Curves are defined as (cluster, index_long) blocks; function grid is index_fn.
# Returns list with: efunctions, evalues, sigma2, argvals, method, obs_prop
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
fgee_fpca_estimate <- function(dx,
                               cname_,
                               index_fn = "yindex.vec",
                               index_long = "time",
                               value_col = "resid",
                               method = c("auto", "face", "sc"),
                               sparse_thresh = 0.70,
                               pve = 0.95,
                               nknots = 15,
                               npc = NULL,
                               var = TRUE,
                               center_sc = FALSE,
                               verbose = TRUE) {

  method <- match.arg(method)

  # COPY ONLY REQUIRED COLUMNS (avoid copying full design matrix)
  needed <- c(cname_, index_fn, index_long, value_col)
  dt <- data.table::as.data.table(dx)[, ..needed]
  dt <- data.table::copy(dt)

  if (!requireNamespace("refund", quietly = TRUE)) {
    stop("Need package 'refund' installed for corr_fn='fpca'.")
  }
  if (!all(c(cname_, index_fn, index_long, value_col) %in% names(dt))) {
    stop("Missing required columns for FPCA: ",
         paste(setdiff(c(cname_, index_fn, index_long, value_col), names(dt)), collapse = ", "))
  }

  # Define curve id = (cluster, index_long)
  dt[, curve_id := interaction(get(cname_), get(index_long), drop = TRUE, lex.order = TRUE)]

  # Full functional grid
  argvals <- sort(unique(dt[[index_fn]]))

  # Map observed points to column positions 1..d (keeps column order stable)
  dt[, idx_pos := match(get(index_fn), argvals)]
  if (anyNA(dt$idx_pos)) stop("FPCA: index_fn values could not be matched to argvals (unexpected).")

  # Wide matrix Y: rows=curves, cols=function grid positions
  # If duplicates exist within (curve_id, idx_pos), we average (robust default).
  wide <- data.table::dcast(
    dt[, .(curve_id, idx_pos, val = get(value_col))],
    curve_id ~ idx_pos,
    value.var = "val",
    fun.aggregate = mean
  )

  Y <- as.matrix(wide[, -1, with = FALSE])

  # Reorder columns numerically
  col_pos <- as.integer(colnames(Y))
  o <- order(col_pos)
  Y <- Y[, o, drop = FALSE]
  col_pos <- col_pos[o]
  argvals_ord <- argvals[col_pos]

  obs_prop <- mean(!is.na(Y))
  any_missing <- anyNA(Y)

  # Choose method
  chosen <- method
  if (method == "auto") {
    # use face only when fully observed (no NA) and not sparse
    if (!any_missing && obs_prop >= sparse_thresh) chosen <- "face" else chosen <- "sc"
  }

  if (verbose) {
    message(sprintf("FPCA: obs_prop=%.3f, missing=%s, method=%s",
                    obs_prop, ifelse(any_missing, "TRUE", "FALSE"), chosen))
  }

  if (chosen == "face") {
    if (any_missing) {
      stop("fpca.face requires dense data (no missing). Use method='sc' or 'auto'.")
    }

    fit <- refund::fpca.face(
      Y = Y,
      argvals = argvals_ord,
      pve = pve,
      npc = npc,
      knots = nknots,
      var = var
    )  # sigma2 is returned when var=TRUE :contentReference[oaicite:4]{index=4}

  } else {
    # fpca.sc supports missing values in an n x d matrix :contentReference[oaicite:5]{index=5}
    fit <- refund::fpca.sc(
      Y = Y,
      argvals = argvals_ord,
      pve = pve,
      npc = npc,
      var = var,
      center = center_sc
    )
  }

  if (is.null(fit$efunctions) || is.null(fit$evalues)) {
    stop("FPCA fit did not return efunctions/evalues as expected.")
  }
  if (is.null(fit$sigma2) || !is.finite(fit$sigma2)) {
    stop("FPCA fit did not return finite sigma2. Ensure var=TRUE.")
  }

  list(
    efunctions = as.matrix(fit$efunctions),
    evalues = as.numeric(fit$evalues),
    sigma2 = as.numeric(fit$sigma2),
    argvals = argvals_ord,
    method = chosen,
    obs_prop = obs_prop
  )
}
