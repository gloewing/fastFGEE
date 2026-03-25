#' @keywords internal
#' @noRd
var.est <- function(di, wi, beta2 = NULL,
                    wi0 = NULL, di0 = NULL, beta0 = NULL,
                    penalty_diag,
                    B = 500,
                    var.type = c("sandwich", "fastboot", "boot"),
                    exact = FALSE,
                    boot.base = c("initial", "final"),
                    return.boot = FALSE,
                    seed = NULL,
                    verbose = FALSE,
                    # backward-compat alias: if callers pass beta=...
                    beta = NULL,
                    # --- NEW: blocking controls for memory-safe fastboot ---
                    block_size = NULL,
                    max_counts_elems = 2e7) {

  var.type  <- match.arg(var.type)
  boot.base <- match.arg(boot.base)

  if (!is.null(seed)) set.seed(seed)

  if (!is.list(di) || !is.list(wi)) stop("di and wi must be lists (one element per cluster).")
  n <- length(wi)
  if (length(di) != n) stop("length(di) must equal length(wi).")
  if (n < 1L) stop("Need at least one cluster.")

  P <- as.matrix(penalty_diag)
  if (!is.numeric(P)) stop("penalty_diag must be numeric.")
  if (nrow(P) != ncol(P)) stop("penalty_diag must be square.")
  p <- nrow(P)

  # Resolve beta2 (allow beta= alias)
  if (is.null(beta2)) beta2 <- beta
  beta2 <- as.numeric(beta2)
  if (length(beta2) != p) stop("beta2 (or beta) must have length = nrow(penalty_diag).")

  # Helper: choose a block size that limits the n x block counts matrix
  .choose_block_size <- function(B, n_rows, block_size, max_counts_elems) {
    if (!is.null(block_size)) {
      bs <- as.integer(block_size)
      if (!is.finite(bs) || bs < 1L) stop("block_size must be a positive integer.")
      return(min(bs, B))
    }
    # default: limit size of counts matrix (n_rows x block)
    bs <- as.integer(floor(max_counts_elems / max(1L, n_rows)))
    bs <- max(1L, min(B, bs))
    bs
  }

  # ---- choose bootstrap base for ONE-STEP bootstrap (exact=FALSE) ----
  # If boot.base="initial" and wi0/di0/beta0 provided, use them.
  if (!isTRUE(exact) && var.type %in% c("boot", "fastboot")) {

    use_initial <- identical(boot.base, "initial") &&
      !is.null(wi0) && !is.null(di0) && !is.null(beta0)

    if (use_initial) {
      wi_boot <- wi0
      di_boot <- di0
      beta_start <- as.numeric(beta0)
      if (length(beta_start) != p) stop("beta0 must have length = nrow(penalty_diag).")
    } else {
      wi_boot <- wi
      di_boot <- di
      beta_start <- beta2
    }

    n0 <- length(wi_boot)
    if (length(di_boot) != n0) stop("length(di_boot) must equal length(wi_boot).")
  }

  # ===========================================================================
  # 1) CENTERED SANDWICH (DEFAULT)
  # ===========================================================================
  if (identical(var.type, "sandwich")) {

    if (n < 2L) stop("var.type='sandwich' requires at least 2 clusters.")

    if (!isTRUE(exact)) {

      # One-step: A = Wbar + P
      Wbar <- Reduce(`+`, wi) / n
      A <- Wbar + P

      Dmat <- do.call(cbind, di)          # p x n

      infl <- solve_pd(A, Dmat)           # p x n
      infl <- infl - rowMeans(infl)

      covhat <- tcrossprod(infl) / (n * (n - 1))
      return(covhat)

    } else {

      # Exact: robust score u_i = d_i - W_i %*% beta2
      Umat <- do.call(
        cbind,
        Map(function(d_i, w_i) d_i - as.numeric(w_i %*% beta2), di, wi)
      ) # p x n

      A_sum <- Reduce(`+`, wi) + P

      infl <- solve_pd(A_sum, Umat)       # p x n
      infl <- infl - rowMeans(infl)

      covhat <- tcrossprod(infl) * (n / (n - 1))
      return(covhat)
    }
  }

  # ===========================================================================
  # 2) FASTBOOT (MEMORY-SAFE BLOCKED + online covariance)
  #    - inverts ONCE and multiplies per block: step = Ainv %*% rhs_block
  # ===========================================================================
  if (identical(var.type, "fastboot")) {

    if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 2L) {
      stop("fastboot requires B >= 2.")
    }
    B <- as.integer(B)

    if (!isTRUE(exact)) {

      n0 <- length(wi_boot)
      if (n0 < 1L) stop("No clusters in bootstrap base.")

      # Fixed A = Wbar + P, invert ONCE
      Wbar <- Reduce(`+`, wi_boot) / n0
      A <- Wbar + P
      Ainv <- solve_pd(A)                 # p x p (ONE inversion)

      # Stack D into p x n0
      Dmat <- do.call(cbind, di_boot)     # p x n0

      # Constant penalty term
      pen_term <- as.numeric(P %*% beta_start) * n0

      # Determine block size
      bs <- .choose_block_size(B, n0, block_size, max_counts_elems)
      if (isTRUE(verbose)) {
        message("fastboot(one-step): n0=", n0, " B=", B, " block_size=", bs)
      }

      prob <- rep.int(1 / n0, n0)

      # --- Two separate code paths to avoid if/else INSIDE the loop ---
      if (isTRUE(return.boot)) {

        boot.beta <- matrix(NA_real_, nrow = B, ncol = p)
        sum_x  <- numeric(p)
        sum_xx <- matrix(0, nrow = p, ncol = p)

        b0 <- 1L
        while (b0 <= B) {
          Bb <- min(bs, B - b0 + 1L)

          counts <- stats::rmultinom(Bb, size = n0, prob = prob)   # n0 x Bb
          Dsum   <- Dmat %*% counts                                # p x Bb

          rhs  <- Dsum - pen_term                                  # p x Bb
          step <- Ainv %*% rhs                                     # p x Bb

          beta_block <- step / n0 + beta_start                      # p x Bb

          # store draws
          boot.beta[b0:(b0 + Bb - 1L), ] <- t(beta_block)

          # online moments for covariance
          sum_x  <- sum_x  + rowSums(beta_block)
          sum_xx <- sum_xx + beta_block %*% t(beta_block)

          b0 <- b0 + Bb
        }

        covhat <- (sum_xx - tcrossprod(sum_x) / B) / (B - 1L)
        return(list(cov = covhat, boot = boot.beta))

      } else {

        sum_x  <- numeric(p)
        sum_xx <- matrix(0, nrow = p, ncol = p)

        b0 <- 1L
        while (b0 <= B) {
          Bb <- min(bs, B - b0 + 1L)

          counts <- stats::rmultinom(Bb, size = n0, prob = prob)   # n0 x Bb
          Dsum   <- Dmat %*% counts                                # p x Bb

          rhs  <- Dsum - pen_term                                  # p x Bb
          step <- Ainv %*% rhs                                     # p x Bb

          beta_block <- step / n0 + beta_start                      # p x Bb

          sum_x  <- sum_x  + rowSums(beta_block)
          sum_xx <- sum_xx + beta_block %*% t(beta_block)

          b0 <- b0 + Bb
        }

        covhat <- (sum_xx - tcrossprod(sum_x) / B) / (B - 1L)
        return(covhat)
      }

    } else {

      # Exact fastboot: keep (sum W + P) fixed; resample only D sums
      Wsum <- Reduce(`+`, wi)
      A <- Wsum + P
      Ainv <- solve_pd(A)                 # p x p (ONE inversion)

      Dmat <- do.call(cbind, di)          # p x n

      bs <- .choose_block_size(B, n, block_size, max_counts_elems)
      if (isTRUE(verbose)) {
        message("fastboot(exact): n=", n, " B=", B, " block_size=", bs)
      }

      prob <- rep.int(1 / n, n)

      if (isTRUE(return.boot)) {

        boot.beta <- matrix(NA_real_, nrow = B, ncol = p)
        sum_x  <- numeric(p)
        sum_xx <- matrix(0, nrow = p, ncol = p)

        b0 <- 1L
        while (b0 <= B) {
          Bb <- min(bs, B - b0 + 1L)

          counts <- stats::rmultinom(Bb, size = n, prob = prob)    # n x Bb
          Dsum   <- Dmat %*% counts                                 # p x Bb

          beta_block <- Ainv %*% Dsum                                # p x Bb

          boot.beta[b0:(b0 + Bb - 1L), ] <- t(beta_block)

          sum_x  <- sum_x  + rowSums(beta_block)
          sum_xx <- sum_xx + beta_block %*% t(beta_block)

          b0 <- b0 + Bb
        }

        covhat <- (sum_xx - tcrossprod(sum_x) / B) / (B - 1L)
        return(list(cov = covhat, boot = boot.beta))

      } else {

        sum_x  <- numeric(p)
        sum_xx <- matrix(0, nrow = p, ncol = p)

        b0 <- 1L
        while (b0 <= B) {
          Bb <- min(bs, B - b0 + 1L)

          counts <- stats::rmultinom(Bb, size = n, prob = prob)    # n x Bb
          Dsum   <- Dmat %*% counts                                 # p x Bb

          beta_block <- Ainv %*% Dsum                                # p x Bb

          sum_x  <- sum_x  + rowSums(beta_block)
          sum_xx <- sum_xx + beta_block %*% t(beta_block)

          b0 <- b0 + Bb
        }

        covhat <- (sum_xx - tcrossprod(sum_x) / B) / (B - 1L)
        return(covhat)
      }
    }
  }

  # ===========================================================================
  # 3) BOOT (TRUE RESAMPLE-AND-RESOLVE)
  # ===========================================================================
  if (identical(var.type, "boot")) {

    if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 2L) {
      stop("boot requires B >= 2.")
    }
    B <- as.integer(B)

    if (!isTRUE(exact)) {

      n0 <- length(wi_boot)
      if (n0 < 1L) stop("No clusters in bootstrap base.")

      pen_term <- as.numeric(P %*% beta_start) * n0
      boot.beta <- matrix(NA_real_, nrow = B, ncol = p)

      for (b in seq_len(B)) {
        samp <- sample.int(n0, n0, replace = TRUE)

        Wbar_b <- Reduce(`+`, wi_boot[samp]) / n0
        Dsum_b <- Reduce(`+`, di_boot[samp])

        step_b <- solve_pd(Wbar_b + P, Dsum_b - pen_term)
        boot.beta[b, ] <- beta_start + as.numeric(step_b) / n0
      }

      covhat <- stats::cov(boot.beta)
      if (return.boot) return(list(cov = covhat, boot = boot.beta))
      return(covhat)

    } else {

      boot.beta <- matrix(NA_real_, nrow = B, ncol = p)

      for (b in seq_len(B)) {
        samp <- sample.int(n, n, replace = TRUE)
        Wsum_b <- Reduce(`+`, wi[samp])
        Dsum_b <- Reduce(`+`, di[samp])
        boot.beta[b, ] <- as.numeric(solve_pd(Wsum_b + P, Dsum_b))
      }

      covhat <- stats::cov(boot.beta)
      if (return.boot) return(list(cov = covhat, boot = boot.beta))
      return(covhat)
    }
  }

  stop("Unreachable var.type branch (this should never happen).")
}

