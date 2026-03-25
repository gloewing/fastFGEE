#' @keywords internal
#' @noRd
.fgee_extract_Ywide_matrix <- function(data, Y_nm) {

  if (!(Y_nm %in% names(data))) stop("Outcome column not found in data: ", Y_nm)

  Y_obj <- data[[Y_nm]]

  # unwrap I(.) / AsIs safely
  if (inherits(Y_obj, "AsIs")) Y_obj <- unclass(Y_obj)

  # 1) already matrix-like
  if (is.matrix(Y_obj) || (is.array(Y_obj) && length(dim(Y_obj)) == 2L)) {
    return(as.matrix(Y_obj))
  }

  # 2) stored as a data.frame/tibble inside a single column (unusual, but happens)
  if (is.data.frame(Y_obj)) {
    return(as.matrix(Y_obj))
  }

  n <- NROW(data)

  # 3) AsIs-wrapped data.frame becomes a *plain list of columns* after unclass()
  #    Detect list-of-columns: length = L (grid), each element length = n (subjects)
  if (is.list(Y_obj) && length(Y_obj) > 1L) {
    len_each <- vapply(Y_obj, length, integer(1))

    if (all(len_each == n) && all(vapply(Y_obj, function(z) is.atomic(z) || is.factor(z), logical(1)))) {
      # coerce each column to numeric and cbind
      cols <- lapply(Y_obj, function(z) as.numeric(z))
      Ymat <- do.call(cbind, cols)
      storage.mode(Ymat) <- "double"
      return(Ymat)
    }
  }

  # 4) list-of-curves: length = n, each entry is a numeric vector (possibly ragged)
  if (is.list(Y_obj) && length(Y_obj) == n) {
    Y_list <- lapply(Y_obj, function(z) {
      if (inherits(z, "AsIs")) z <- unclass(z)
      if (is.matrix(z) || (is.array(z) && length(dim(z)) == 2L)) z <- as.numeric(z)
      as.numeric(z)
    })
    if (all(vapply(Y_list, is.numeric, logical(1)))) {
      L <- max(vapply(Y_list, length, integer(1)))
      Ymat <- matrix(NA_real_, nrow = n, ncol = L)
      for (i in seq_len(n)) {
        yi <- Y_list[[i]]
        if (length(yi)) Ymat[i, seq_along(yi)] <- yi
      }
      return(Ymat)
    }
  }

  # 5) scalar outcome vector => treat as 1-col matrix
  if (is.atomic(Y_obj) && length(Y_obj) == n) {
    return(matrix(as.numeric(Y_obj), ncol = 1))
  }

  stop(
    "Could not coerce data[['", Y_nm, "']] into a numeric matrix.\n",
    "Class(data[['", Y_nm, "']]) = ", paste(class(data[[Y_nm]]), collapse = "/"), "\n",
    "Try normalizing your wide response before calling fgee/pffr:\n",
    "  data[['", Y_nm, "']] <- I(as.matrix(data[['", Y_nm, "']]))\n"
  )
}

#' @keywords internal
#' @noRd
.fgee_build_long_data <- function(data,
                                  fit_pffr,
                                  formula,
                                  cluster,
                                  time = NULL,
                                  index_fn = "yindex.vec",
                                  index_long = "time",
                                  y_col = "Y",
                                  cluster_col = "cluster",
                                  cname_col = "cname_",
                                  canonical_order = FALSE,
                                  order_cols = NULL,
                                  prefer_fit_X = TRUE,
                                  check_alignment = TRUE) {


  if (!inherits(formula, "formula")) stop("`formula` must be a formula.")
  if (is.null(fit_pffr)) stop("fit_pffr cannot be NULL.")
  if (!is.character(cluster) || length(cluster) != 1L) stop("cluster must be a single column name.")
  if (!is.data.frame(data)) stop("data must be a data.frame (wide). data.table is allowed but will be treated as a list-like object.")

  # IMPORTANT: do NOT convert wide data to data.table here.
  nm <- names(data)
  if (!(cluster %in% nm)) stop("cluster variable not found in data: ", cluster)

  # outcome name (wide functional response matrix stored as a column)
  Y_nm <- all.vars(formula)[1]
  if (!(Y_nm %in% nm)) stop("Outcome variable not found in data: ", Y_nm)

  # Pull wide response safely (avoid printing it)
  Y_wide <- .fgee_extract_Ywide_matrix(data, Y_nm)

  # observed lengths per subject (matches your old Y_len logic)
  Y_len <- rowSums(!is.na(Y_wide))
  ID <- rep.int(data[[cluster]], times = Y_len)

  # Allow a few common storage patterns
  if (is.data.frame(Y_wide)) Y_wide <- as.matrix(Y_wide)

  # If a tibble/list-column of curves was used, convert to matrix (ragged supported)
  if (is.list(Y_wide) && length(Y_wide) == nrow(data) && all(vapply(Y_wide, is.numeric, logical(1)))) {
    L <- max(vapply(Y_wide, length, integer(1)))
    Ymat <- matrix(NA_real_, nrow = nrow(data), ncol = L)
    for (i in seq_len(nrow(data))) {
      yi <- Y_wide[[i]]
      Ymat[i, seq_along(yi)] <- yi
    }
    Y_wide <- Ymat
  }

  if (!is.matrix(Y_wide)) {
    stop("Expected data[['", Y_nm, "']] to be a matrix (or AsIs-wrapped matrix). Got: ", paste(class(data[[Y_nm]]), collapse = "/"))
  }

  n <- nrow(Y_wide)
  # observed lengths per subject (matches your old Y_len)
  Y_len <- rowSums(!is.na(Y_wide))
  if (any(!is.finite(Y_len))) stop("Non-finite Y_len encountered.")

  # ID and time vectors repeated by Y_len (matches your old logic)
  ID <- rep.int(data[[cluster]], times = Y_len)

  if (is.null(time)) {
    time_vec <- rep.int(1, length(ID))
  } else if (is.character(time) && length(time) == 1L) {
    if (!(time %in% nm)) stop("time column not found in data: ", time)
    time_vec <- rep.int(data[[time]], times = Y_len)
  } else {
    if (length(time) != nrow(data)) stop("time must be NULL, a column name, or a vector of length nrow(data).")
    time_vec <- rep.int(time, times = Y_len)
  }

  # ---- Pull long response + index from the fitted object (preferred for alignment) ----
  mf <- fit_pffr$model
  if (is.null(mf)) stop("fit_pffr$model is NULL; cannot build aligned long data.")

  # response column inside fit$model
  y_fit <- NULL
  if (Y_nm %in% names(mf)) {
    y_fit <- mf[[Y_nm]]
  } else if ("Y" %in% names(mf)) {
    y_fit <- mf[["Y"]]
  } else {
    stop("Could not find response column in fit_pffr$model (looked for '", Y_nm, "' and 'Y').")
  }
  y_fit <- as.numeric(y_fit)

  # functional index inside fit$model
  if (index_fn %in% names(mf)) {
    idx_fn_fit <- mf[[index_fn]]
  } else if ("yindex.vec" %in% names(mf)) {
    idx_fn_fit <- mf[["yindex.vec"]]
  } else {
    stop("Could not find index_fn='", index_fn, "' in fit_pffr$model.")
  }
  idx_fn_fit <- as.numeric(idx_fn_fit)

  # design matrix (do NOT rebuild formulas from colnames)
  MM <- NULL
  if (isTRUE(prefer_fit_X)) {
    MM <- tryCatch(fit_pffr$X, error = function(e) NULL)
  }
  if (is.null(MM)) {
    MM <- suppressWarnings(stats::model.matrix(fit_pffr))
  }
  if (inherits(MM, "Matrix")) MM <- as.matrix(MM)
  if (!is.matrix(MM)) stop("Could not obtain a matrix design for fit_pffr (X/model.matrix).")

  # ---- Alignment checks (same spirit as your old guard) ----
  if (length(ID) != length(y_fit)) {
    stop(
      "Some rows/curves were dropped in pffr() due to NA covariates or exclusions.\n",
      "Constructed long ID length = ", length(ID),
      " but fit_pffr long response length = ", length(y_fit), ".\n",
      "Remove/handle those wide rows before fitting, or build ID/time from fit_pffr$model instead."
    )
  }
  if (length(idx_fn_fit) != length(y_fit)) stop("Length mismatch: fit yindex length != fit y length.")
  if (nrow(MM) != length(y_fit)) stop("Row mismatch: nrow(design) != length(y_fit).")

  if (isTRUE(check_alignment)) {
    y2 <- tryCatch(fit_pffr$y, error = function(e) NULL)
    if (!is.null(y2) && length(y2) == length(y_fit)) {
      ii <- unique(pmin(length(y_fit), c(1:5, round(seq(1, length(y_fit), length.out = 10)))))
      if (!isTRUE(all.equal(as.numeric(y2[ii]), y_fit[ii], tolerance = 1e-10))) {
        stop("Alignment check failed: fit_pffr$y and fit_pffr$model response differ at sampled indices.")
      }
    }
  }

  # ---- Build LONG data.table (safe: atomic columns only) ----
  dx <- data.table::data.table()
  dx[, (y_col) := y_fit]
  dx[, (cluster_col) := ID]
  dx[, (index_fn) := idx_fn_fit]
  dx[, (index_long) := time_vec]
  dx <- cbind(dx, data.table::as.data.table(MM))
  dx[, (cname_col) := get(cluster_col)]

  if (isTRUE(canonical_order)) {
    if (is.null(order_cols)) order_cols <- c(cluster_col, index_long, index_fn)
    data.table::setorderv(dx, order_cols)
  }

  list(
    dx = dx,
    X_ = colnames(MM),
    namesd = colnames(MM),
    yindex.vec = sort(unique(idx_fn_fit))
  )
}
