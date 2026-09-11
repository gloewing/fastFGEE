# Exact fast cluster cross-validation workspaces and optimizers ----------------

# Canonical family key for the fastK path.
#
# A fitted mgcv extended family carries its estimated parameter in the family
# string, so `mgcv::nb()` becomes "Negative Binomial(2.338)" and `mgcv::betar()`
# becomes "Beta regression(13.227)".  The digit-preserving normalisation below
# would turn those into "negativebinomial2338" and "betaregression13227", which
# no fixed key can match -- and `.fgee_fastk_family_code()` treats an unmatched
# key as a silent fall-through to the R path, so the failure would be invisible.
#
# Negative binomial and beta are therefore delegated to
# `.fgee_nuisance_family_key()`, the one normaliser in the package that
# lowercases before stripping and matches on patterns.  Its `gamma` answer is
# deliberately *not* used: that regex is `grepl("gamma", z)`, which would
# collapse "quasigamma" into "gamma" and erase a distinction the fastK gate and
# the kernel family codes rely on.
.fgee_family_key <- function(family) {
  fam <- if (inherits(family, "family")) family$family else as.character(family)[1L]
  key <- .fgee_nuisance_family_key(family)
  if (key %in% c("negbinomial", "beta")) return(key)
  tolower(gsub("[^A-Za-z0-9]", "", fam))
}

.fgee_initial_lambda <- function(fit.initial, penalty_setup) {
  sp <- fit.initial$sp
  if (is.list(sp)) sp <- sp[[1L]]
  sp <- as.numeric(sp)
  if (!length(sp)) sp <- rep(1, max(1L, penalty_setup$n_sm))
  sp[!is.finite(sp) | sp <= 0] <- 1
  sp
}

.fgee_make_cluster_folds <- function(ids, K = 10L, seed = 1L, sets = NULL) {
  ids <- as.character(ids)
  n <- length(ids)
  K <- min(as.integer(K), n)
  if (K < 2L) stop("At least two folds are required.")

  if (is.null(sets)) {
    set.seed(seed)
    indices <- seq_len(n)
    sets <- split(indices[sample.int(n)], sort(rank(indices) %% K))
  } else {
    sets <- lapply(sets, as.integer)
    K <- length(sets)
  }

  used <- sort(unlist(sets, use.names = FALSE))
  if (!identical(used, seq_len(n))) {
    stop("sets must partition all cluster indices exactly once.")
  }

  fold_id <- integer(n)
  for (k in seq_len(K)) fold_id[sets[[k]]] <- k
  names(fold_id) <- ids

  list(sets = sets, fold_id = fold_id, K = K, ids = ids, seed = seed)
}

.fgee_fastk_gaussian_stats <- function(dx, namesd, working, folds,
                                        cname_ = "cname_") {
  p <- length(namesd)
  K <- folds$K
  N <- working$N

  if (!is.null(working$gaussian_cv)) {
    old <- working$gaussian_cv$fold_id[working$cluster_id]
    new <- folds$fold_id[working$cluster_id]
    if (identical(as.integer(old), as.integer(new))) return(working$gaussian_cv)
  }

  C <- array(0, c(p, p, K))
  cc <- matrix(0, p, K)
  ss <- numeric(K)
  cid <- as.character(dx[[cname_]])
  rows <- split(seq_len(nrow(dx)), cid)
  rows <- rows[working$cluster_id]

  for (i in seq_len(N)) {
    ii <- rows[[i]]
    X <- as.matrix(dx[ii, ..namesd])
    y <- as.numeric(dx$Y[ii])
    sc <- 1 / (N * length(ii))
    k <- folds$fold_id[i]
    C[, , k] <- C[, , k] + sc * crossprod(X)
    cc[, k] <- cc[, k] + sc * as.numeric(crossprod(X, y))
    ss[k] <- ss[k] + sc * sum(y^2)
  }

  list(
    C = C,
    c = cc,
    s = ss,
    fold_id = setNames(folds$fold_id, working$cluster_id),
    fold_cluster_count = vapply(folds$sets, length, integer(1))
  )
}

#' Prepare the exact fastK objective
#' @keywords internal
#' @noRd
fgee_fastk_prepare <- function(working,
                               dx,
                               namesd,
                               fit.initial,
                               K = 10L,
                               seed = 1L,
                               sets = NULL,
                               exact = FALSE,
                               memory = c("balanced", "speed", "lowmem"),
                               clip_prob = 1e-6,
                               clip_y = 1e-6,
                               nuisance = NULL,
                               cname_ = "cname_") {
  if (!inherits(working, "fgee_working_stats")) {
    stop("working must be an fgee_working_stats object.")
  }
  memory <- match.arg(memory)
  D <- .fgee_working_dmat(working, exact = isTRUE(exact), required = TRUE)
  folds <- .fgee_make_cluster_folds(working$cluster_id, K = K, seed = seed,
                                    sets = sets)

  beta0 <- as.numeric(fit.initial$coefficients)
  p <- length(beta0)
  beta_mat <- matrix(beta0, p, folds$K)
  total_rows <- sum(working$cluster_size)

  d_hold <- vapply(
    folds$sets,
    function(ii) rowSums(D[, ii, drop = FALSE]),
    numeric(p)
  )
  if (folds$K == 1L) d_hold <- matrix(d_hold, p, 1L)
  holdout_rows <- vapply(
    folds$sets,
    function(ii) sum(working$cluster_size[ii]),
    numeric(1)
  )

  d_total <- if (isTRUE(exact)) working$d_exact_sum else working$d_sum
  if (is.null(d_total) || length(d_total) != p) {
    stop("Required cluster score totals are unavailable for fastK preparation.")
  }
  D_train <- matrix(d_total, p, folds$K) - d_hold
  D_train <- sweep(D_train, 2L, total_rows / (total_rows - holdout_rows), "*")

  ps <- penalty_setup(fit.initial, unpenalized = fit.initial$nsdf)
  frem <- .fgee_initial_lambda(fit.initial, ps)
  q <- length(frem)
  components <- penalty_components_from_setup(ps, q)

  key <- .fgee_family_key(fit.initial$family)
  link <- tolower(fit.initial$family$link)
  gaussian <- key %in% c("gaussian", "normal") && link == "identity"

  supported <- gaussian ||
    (key %in% c("binomial", "quasibinomial") && link == "logit") ||
    (key %in% c("poisson", "quasipoisson") && link == "log") ||
    (key %in% c("gamma", "quasigamma") && link == "log") ||
    (identical(key, "negbinomial") && link == "log") ||
    (identical(key, "beta") && link == "logit")
  if (!supported) {
    stop(
      "Analytic fastK supports Gaussian/identity, Binomial/logit, ",
      "Poisson/log, Gamma/log, Negative Binomial/log and Beta/logit. Got ",
      fit.initial$family$family, "/", fit.initial$family$link, "."
    )
  }

  # Beta regression is defined on the open unit interval. Values at exactly 0 or
  # 1 are clipped by the loss (clip_y); values outside [0, 1] are a data error,
  # and catching them here is far clearer than an NaN surfacing from log().
  if (identical(key, "beta")) {
    yv <- data.table::as.data.table(dx)[["Y"]]
    if (is.null(yv)) stop("Beta regression fastK needs a 'Y' column.",
                          call. = FALSE)
    rg <- suppressWarnings(range(as.numeric(yv)))
    if (!all(is.finite(rg)) || rg[1L] < 0 || rg[2L] > 1) {
      stop("Beta regression requires responses in [0, 1]; observed range [",
           format(rg[1L]), ", ", format(rg[2L]), "].", call. = FALSE)
    }
  }

  # Scalar nuisance for the criterion.
  #
  # This must be the value that produced `working`, not a fresh estimate.
  # Bhat[, k] in fgee_fastk_score_grad() is built from prep$Wbar and prep$D,
  # both functions of the working variance -- mu + mu^2/theta for negative
  # binomial, mu(1-mu)/(1+phi) for beta -- so evaluating those coefficients
  # under a *different* nuisance would stop being a held-out predictive loss for
  # the model actually being fitted, and the resulting bias would follow neither
  # value. So read it from working$nuisance and never estimate it here.
  #
  # `nuisance` overrides that, for the sensitivity study only; it deliberately
  # mirrors `clip_prob`, an argument the dispatch layer never sets.
  nu_key <- .fgee_nuisance_family_key(fit.initial$family)
  nu_param <- .fgee_nuisance_parameter_name(nu_key)
  nu_value <- NA_real_
  nu_source <- NA_character_
  if (!is.null(nu_param)) {
    if (!is.null(nuisance)) {
      nu_value <- as.numeric(nuisance)[1L]
      nu_source <- "override"
    } else {
      v <- .fgee_nuisance_value(working, nu_param, NULL)
      if (.fgee_positive_scalar(v)) {
        nu_value <- as.numeric(v)[1L]
        nu_source <- "working$nuisance"
      }
    }
    # Gamma's dispersion cancels out of the criterion (it enters the
    # log-likelihood multiplicatively, plus eta-free terms), so a missing value
    # is harmless there. For negative binomial and beta it does not cancel, and
    # a missing value is a precondition violation rather than something to
    # silently default.
    if (nu_key %in% c("negbinomial", "beta") && !.fgee_positive_scalar(nu_value)) {
      stop(
        "Analytic fastK for ", fit.initial$family$family, "/",
        fit.initial$family$link, " requires a positive ", nu_param,
        "; none was available from the working statistics. Supply it through ",
        "the family (see ?mgcv::nb, ?mgcv::betar) or pass nuisance=.",
        call. = FALSE
      )
    }
  }

  out <- list(
    working = working,
    fit.initial = fit.initial,
    ncl = working$N,
    p = p,
    q = q,
    K = folds$K,
    folds = folds,
    beta0 = beta0,
    beta_mat = beta_mat,
    D = D_train,
    Wbar = working$W_bar,
    ps = ps,
    components = components$S,
    lambda_frem = frem,
    family_key = key,
    link = link,
    gaussian = gaussian,
    clip_prob = clip_prob,
    clip_y = clip_y,
    nuisance_value = nu_value,
    nuisance_source = nu_source,
    memory = memory,
    cname_ = cname_,
    namesd = namesd,
    exact = isTRUE(exact)
  )

  if (gaussian) {
    out$gaussian_stats <- .fgee_fastk_gaussian_stats(
      dx, namesd, working, folds, cname_ = cname_
    )
  } else {
    dd <- data.table::as.data.table(dx)
    cid <- as.character(dd[[cname_]])
    row_map <- split(seq_len(nrow(dd)), cid)
    row_map <- row_map[working$cluster_id]

    if (memory == "lowmem") {
      out$dx <- dd
      out$row_map <- row_map
    } else if (memory == "speed") {
      XX <- YY <- WW <- vector("list", folds$K)
      for (k in seq_len(folds$K)) {
        rr <- unlist(row_map[folds$sets[[k]]], use.names = FALSE)
        Xk <- as.matrix(dd[rr, ..namesd])
        yk <- as.numeric(dd$Y[rr])
        ids_k <- working$cluster_id[folds$sets[[k]]]
        gk <- match(cid[rr], ids_k)
        nk <- working$cluster_size[folds$sets[[k]]]
        XX[[k]] <- Xk
        YY[[k]] <- yk
        WW[[k]] <- 1 / (working$N * nk[gk])
      }
      out$XX <- XX
      out$YY <- YY
      out$row_weight <- WW
    } else {
      row_blocks <- vector("list", folds$K)
      weight_blocks <- vector("list", folds$K)
      block_length <- integer(folds$K)
      for (k in seq_len(folds$K)) {
        ids_k <- working$cluster_id[folds$sets[[k]]]
        rr <- unlist(row_map[ids_k], use.names = FALSE)
        gk <- match(cid[rr], ids_k)
        nk <- working$cluster_size[folds$sets[[k]]]
        row_blocks[[k]] <- rr
        weight_blocks[[k]] <- 1 / (working$N * nk[gk])
        block_length[k] <- length(rr)
      }
      row_order <- unlist(row_blocks, use.names = FALSE)
      out$X_eval <- as.matrix(dd[row_order, ..namesd])
      out$Y_eval <- as.numeric(dd$Y[row_order])
      out$W_eval <- unlist(weight_blocks, use.names = FALSE)
      out$fold_ends <- cumsum(block_length)
      out$fold_starts <- c(1L, head(out$fold_ends, -1L) + 1L)
    }
  }

  class(out) <- c("fgee_fastk_workspace", "list")
  out
}

.fgee_fastk_fold_data <- function(prep, k) {
  if (!is.null(prep$XX)) {
    return(list(X = prep$XX[[k]], y = prep$YY[[k]], w = prep$row_weight[[k]]))
  }
  if (!is.null(prep$X_eval)) {
    ii <- prep$fold_starts[k]:prep$fold_ends[k]
    return(list(
      X = prep$X_eval[ii, , drop = FALSE],
      y = prep$Y_eval[ii],
      w = prep$W_eval[ii]
    ))
  }

  ids <- prep$working$cluster_id[prep$folds$sets[[k]]]
  rr <- unlist(prep$row_map[ids], use.names = FALSE)
  cid <- as.character(prep$dx[[prep$cname_]][rr])
  g <- match(cid, ids)
  nk <- prep$working$cluster_size[prep$folds$sets[[k]]]
  list(
    X = as.matrix(prep$dx[rr, prep$namesd, with = FALSE]),
    y = as.numeric(prep$dx$Y[rr]),
    w = 1 / (prep$working$N * nk[g])
  )
}

.fgee_fastk_lowmem_fold <- function(prep, k, beta, need_gradient = TRUE) {
  ids_idx <- prep$folds$sets[[k]]
  loss <- 0
  grad <- if (isTRUE(need_gradient)) numeric(prep$p) else NULL

  for (i in ids_idx) {
    id <- prep$working$cluster_id[i]
    rr <- prep$row_map[[id]]
    if (is.null(rr) || !length(rr)) stop("Missing rows for cluster '", id, "'.")

    eta <- numeric(length(rr))
    for (j in seq_len(prep$p)) {
      bj <- beta[j]
      if (is.finite(bj) && bj != 0) {
        eta <- eta + prep$dx[[prep$namesd[j]]][rr] * bj
      }
    }

    y <- as.numeric(prep$dx$Y[rr])
    ld <- .fgee_fastk_loss_deta(prep, y, eta)
    w <- 1 / (prep$ncl * prep$working$cluster_size[i])
    loss <- loss + w * sum(ld$loss)

    if (isTRUE(need_gradient)) {
      r <- w * ld$deta
      for (j in seq_len(prep$p)) {
        grad[j] <- grad[j] + sum(prep$dx[[prep$namesd[j]]][rr] * r)
      }
    }
  }

  list(loss = loss, gradient_beta = grad)
}

.fgee_fastk_loss_deta <- function(prep, y, eta) {
  key <- prep$family_key

  if (key %in% c("gaussian", "normal")) {
    e <- eta - y
    return(list(loss = e^2, deta = 2 * e))
  }
  if (key %in% c("binomial", "quasibinomial")) {
    p <- stats::plogis(eta)
    eps <- prep$clip_prob
    pc <- pmin(pmax(p, eps), 1 - eps)
    loss <- -(y * log(pc) + (1 - y) * log(1 - pc))
    active <- p > eps & p < 1 - eps
    deta <- numeric(length(y))
    deta[active] <- p[active] - y[active]
    return(list(loss = loss, deta = deta))
  }
  if (key %in% c("poisson", "quasipoisson")) {
    mu <- exp(eta)
    mc <- pmax(mu, 1e-12)
    loss <- mc - y * log(mc)
    deta <- numeric(length(y))
    active <- mu > 1e-12
    deta[active] <- mu[active] - y[active]
    return(list(loss = loss, deta = deta))
  }
  if (key %in% c("gamma", "quasigamma")) {
    mu <- exp(eta)
    mc <- pmax(mu, 1e-12)
    loss <- log(mc) + y / mc - 1
    deta <- numeric(length(y))
    active <- mu > 1e-12
    deta[active] <- 1 - y[active] / mu[active]
    return(list(loss = loss, deta = deta))
  }
  if (identical(key, "negbinomial")) {
    # NB2: Var(Y) = mu + mu^2/theta.  Dropping the eta-free terms of the
    # log-likelihood leaves (theta + y) log(theta + mu) - y log(mu), and with a
    # log link
    #     dloss/deta = mu (theta + y)/(theta + mu) - y = theta (mu - y)/(mu + theta).
    #
    # The expression shape deliberately matches the legacy path at R/cv.R:346,
    # which lets the staged-versus-legacy agreement test cross-check the score
    # itself against an independently written implementation.  It keeps the
    # eta-free (theta + y) log(theta) term, so the score is not directly
    # comparable to the Poisson score as theta grows; the large-theta limit is
    # therefore asserted on the *gradient*, where that term does not appear and
    # theta (mu - y)/(mu + theta) -> mu - y exactly.
    #
    # Clipping follows the sibling log-link families (1e-12 floor, and the
    # zero-derivative-when-clipped convention) rather than cv.R's 1e-10.  No
    # upper clip: dloss/deta is bounded on [-y, theta] for every mu > 0 and the
    # loss is finite for any finite mu, so an upper clip would only add a
    # plateau on which the gradient vanishes at the loss maximum.
    th <- prep$nuisance_value
    if (!.fgee_positive_scalar(th)) {
      stop("The negative-binomial fastK loss requires a positive theta.",
           call. = FALSE)
    }
    mu <- exp(eta)
    mc <- pmax(mu, 1e-12)
    loss <- (th + y) * log(th + mc) - y * log(mc)
    deta <- numeric(length(y))
    active <- mu > 1e-12
    deta[active] <- th * (mu[active] - y[active]) / (mu[active] + th)
    return(list(loss = loss, deta = deta))
  }
  if (identical(key, "beta")) {
    # Beta regression with precision phi: Var(Y) = mu (1 - mu)/(1 + phi), and
    # Y ~ Beta(a, b) with a = mu phi, b = (1 - mu) phi.  With a logit link
    #     dloss/deta = phi mu (1 - mu) [digamma(a) - digamma(b) - log y + log(1 - y)].
    #
    # The expression below is the one at R/cv.R:361, term for term, and that is
    # load-bearing rather than stylistic.  Individual lgamma terms reach ~180
    # while the loss itself passes through zero, so two algebraically identical
    # arrangements of this sum differ by ~2e-10 -- two hundred times the
    # tolerance the compiled-kernel agreement tests use.  In particular do not
    # simplify lgamma(a) + lgamma(b) - lgamma(a + b) using a + b == phi, and use
    # log(1 - yb) rather than log1p(-yb).  Any future compiled mirror of this
    # branch must copy the arrangement, not just the constants.
    ph <- prep$nuisance_value
    if (!.fgee_positive_scalar(ph)) {
      stop("The beta-regression fastK loss requires a positive precision.",
           call. = FALSE)
    }
    eps <- prep$clip_prob
    ey <- .fgee_or(prep$clip_y, 1e-6)
    mu <- stats::plogis(eta)
    mc <- pmin(pmax(mu, eps), 1 - eps)
    yb <- pmin(pmax(y, ey), 1 - ey)
    a <- mc * ph
    b <- (1 - mc) * ph
    lyb <- log(yb)
    l1yb <- log(1 - yb)
    loss <- lgamma(a) + lgamma(b) - lgamma(a + b) -
      (a - 1) * lyb - (b - 1) * l1yb
    deta <- numeric(length(y))
    active <- mu > eps & mu < 1 - eps
    deta[active] <- ph * mu[active] * (1 - mu[active]) *
      (digamma(a[active]) - digamma(b[active]) - lyb[active] + l1yb[active])
    return(list(loss = loss, deta = deta))
  }
  stop("Unsupported family in fastK loss.")
}

#' Evaluate exact fastK and its analytic gradient
#' @keywords internal
#' @noRd
fgee_fastk_score_grad <- function(prep, lambda, need_gradient = TRUE,
                                  return_fold = FALSE) {
  lambda <- as.numeric(lambda)
  if (length(lambda) != prep$q || any(!is.finite(lambda)) || any(lambda <= 0)) {
    return(list(score = Inf, gradient = rep(NA_real_, prep$q)))
  }

  P <- penalty_from_setup(prep$ps, lambda)
  if (isTRUE(prep$exact)) {
    fac <- .fgee_factor_pd(prep$Wbar + P / prep$ncl)
    Bhat <- fac$solve(prep$D / prep$ncl)
  } else {
    fac <- .fgee_factor_pd(prep$Wbar + P)
    rhs <- prep$D - (P %*% prep$beta_mat) * prep$ncl
    Bhat <- prep$beta_mat + fac$solve(rhs) / prep$ncl
  }

  fold_contribution <- numeric(prep$K)
  Gbeta <- if (need_gradient) matrix(0, prep$p, prep$K) else NULL

  if (isTRUE(prep$gaussian)) {
    gs <- prep$gaussian_stats
    for (k in seq_len(prep$K)) {
      b <- Bhat[, k]
      Ck <- gs$C[, , k]
      ck <- gs$c[, k]
      fold_contribution[k] <- gs$s[k] - 2 * sum(ck * b) +
        sum(b * as.numeric(Ck %*% b))
      if (need_gradient) Gbeta[, k] <- 2 * (as.numeric(Ck %*% b) - ck)
    }
  } else if (identical(prep$memory, "lowmem")) {
    for (k in seq_len(prep$K)) {
      ev_k <- .fgee_fastk_lowmem_fold(
        prep, k, beta = Bhat[, k], need_gradient = need_gradient
      )
      fold_contribution[k] <- ev_k$loss
      if (need_gradient) Gbeta[, k] <- ev_k$gradient_beta
    }
  } else {
    for (k in seq_len(prep$K)) {
      fd <- .fgee_fastk_fold_data(prep, k)
      eta <- as.numeric(fd$X %*% Bhat[, k])
      ld <- .fgee_fastk_loss_deta(prep, fd$y, eta)
      fold_contribution[k] <- sum(fd$w * ld$loss)
      if (need_gradient) {
        Gbeta[, k] <- as.numeric(crossprod(fd$X, fd$w * ld$deta))
      }
    }
  }

  score <- sum(fold_contribution)
  fold_weight <- vapply(prep$folds$sets, length, integer(1)) / prep$ncl
  fold_mean <- fold_contribution / fold_weight

  if (!need_gradient) {
    ans <- list(score = score, beta_folds = Bhat)
    if (return_fold) ans$fold_score <- fold_mean
    return(ans)
  }

  grad <- numeric(prep$q)
  for (j in seq_len(prep$q)) {
    dP <- log(10) * lambda[j] * prep$components[[j]]
    if (isTRUE(prep$exact)) dP <- dP / prep$ncl
    dB <- -fac$solve(dP %*% Bhat)
    grad[j] <- sum(Gbeta * dB)
  }

  ans <- list(score = score, gradient = grad, beta_folds = Bhat, Gbeta = Gbeta)
  if (return_fold) ans$fold_score <- fold_mean
  ans
}

.fgee_refine_grid <- function(lambda, multipliers, zero_floor = 1e-8) {
  lambda <- pmax(as.numeric(lambda), zero_floor)
  q <- length(lambda)
  if (is.list(multipliers)) {
    if (length(multipliers) == 1L) multipliers <- rep(multipliers, q)
    if (length(multipliers) != q) stop("Multiplier list has the wrong length.")
    vals <- lapply(seq_len(q), function(j) lambda[j] * as.numeric(multipliers[[j]]))
  } else if (is.matrix(multipliers) || is.data.frame(multipliers)) {
    multipliers <- as.data.frame(multipliers)
    if (ncol(multipliers) == 1L) multipliers <- multipliers[rep(1L, q)]
    if (ncol(multipliers) != q) stop("Multiplier matrix has the wrong number of columns.")
    vals <- lapply(seq_len(q), function(j) lambda[j] * as.numeric(multipliers[[j]]))
  } else {
    vals <- lapply(lambda, function(x) x * as.numeric(multipliers))
  }
  as.matrix(expand.grid(vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
}

.fgee_prepare_staged_grids <- function(cv.grid, lambda_frem) {
  if (!is.list(cv.grid)) cv.grid <- list(cv.grid)
  if (!length(cv.grid)) stop("cv.grid must contain at least one stage.")
  out <- cv.grid
  out[[1L]] <- do.call(cbind, lapply(lambda_frem, function(x) x * cv.grid[[1L]]))
  out
}

#' Staged exact fastK using the compact workspace
#' @keywords internal
#' @noRd
fgee_tune_fastk_staged <- function(prep, cv.grid, zero_floor = 1e-8,
                                   verbose = TRUE) {
  grids <- .fgee_prepare_staged_grids(cv.grid, prep$lambda_frem)
  scores <- fold_scores <- stage_se <- vector("list", length(grids))
  final_grid <- NULL

  for (stage in seq_along(grids)) {
    if (verbose) message("fast k-fold stage ", stage, " of ", length(grids))
    grid <- as.matrix(grids[[stage]])
    sc <- numeric(nrow(grid))
    fs <- matrix(NA_real_, nrow(grid), prep$K)

    for (g in seq_len(nrow(grid))) {
      ev <- fgee_fastk_score_grad(prep, grid[g, ], need_gradient = FALSE,
                                  return_fold = TRUE)
      sc[g] <- ev$score
      fs[g, ] <- ev$fold_score
    }

    scores[[stage]] <- sc
    fold_scores[[stage]] <- fs
    stage_se[[stage]] <- apply(fs, 1L, stats::sd) / sqrt(prep$K)
    best <- which.min(sc)
    if (stage < length(grids)) {
      grids[[stage + 1L]] <- .fgee_refine_grid(
        grid[best, ], grids[[stage + 1L]], zero_floor = zero_floor
      )
    } else {
      final_grid <- grid
      final_score <- sc
      final_fold <- fs
    }
  }

  # Match the published/legacy minimum rule: when several final candidates
  # have exactly the same loss, prefer the smoother candidate (largest total
  # penalty trace).  The tie rule matters on the very flat fastK surfaces that
  # motivated the continuous optimizer.
  i_min <- which.min(final_score)
  tied <- which(final_score == final_score[i_min])
  if (length(tied) > 1L) {
    trace_pen <- vapply(tied, function(i) {
      sum(diag(penalty_from_setup(prep$ps, final_grid[i, ])))
    }, numeric(1))
    i <- tied[which.max(trace_pen)]
  } else {
    i <- i_min
  }
  lambda <- as.numeric(final_grid[i, ])
  list(
    method = "fastk_staged",
    lambda = lambda,
    lambda.star = matrix(lambda, nrow = 1L),
    score = final_score[i],
    mse = final_score,
    se = stage_se[[length(stage_se)]],
    mse.ls = scores,
    se.ls = stage_se,
    fold_scores = fold_scores,
    grid = final_grid,
    argmin = i_min,
    argfinal = i,
    lambda.star.min = matrix(as.numeric(final_grid[i_min, ]), nrow = 1L),
    penalty_mat = penalty_from_setup(prep$ps, lambda),
    evaluations = sum(vapply(scores, length, integer(1))),
    convergence = 0L,
    boundary = FALSE,
    workspace = prep
  )
}

.fgee_default_ray <- function() {
  sort(unique(c(
    10^(-6:4),
    10^(-6:4) * 0.25,
    10^(-6:4) * 0.50,
    10^(-6:4) * 0.75,
    10^(-6:4) * 0.90
  )))
}

.fgee_fastk_start_grid <- function(prep,
                                    qreml_lambda = NULL,
                                    strategy = c("qreml", "robust", "compact"),
                                    bound = 8) {
  strategy <- match.arg(strategy)
  q <- prep$q
  base <- prep$lambda_frem
  base[!is.finite(base) | base <= 0] <- 1

  starts <- matrix(0, nrow = 1L, ncol = q)
  ray <- log10(.fgee_default_ray())
  ray <- ray[abs(ray) <= bound]
  starts <- rbind(starts, do.call(rbind, lapply(ray, function(a) rep(a, q))))

  if (!is.null(qreml_lambda)) {
    uq <- log10(as.numeric(qreml_lambda) / base)
    uq <- pmin(pmax(uq, -bound), bound)
    if (all(is.finite(uq))) starts <- rbind(starts, uq)
  }

  axis <- if (strategy == "robust") c(-4, -2, 2, 4) else c(-2, 2)
  if (strategy != "compact") {
    for (j in seq_len(q)) {
      for (a in axis) {
        u <- rep(0, q)
        u[j] <- a
        starts <- rbind(starts, u)
      }
    }
  }

  if (strategy == "robust" && q <= 5L) {
    starts <- rbind(starts, as.matrix(expand.grid(rep(list(c(-3, 3)), q))))
  }

  starts <- unique(starts)
  starts[apply(abs(starts) <= bound, 1L, all), , drop = FALSE]
}

#' Continuous analytic-gradient exact fastK
#' @keywords internal
#' @noRd
fgee_tune_fastk_grad <- function(prep,
                                 qreml_lambda = NULL,
                                 start_strategy = c("qreml", "robust", "compact"),
                                 bound = 8,
                                 n_starts = 4L,
                                 maxit = 100L,
                                 factr = 1e8,
                                 pgtol = 1e-7,
                                 verbose = FALSE) {
  start_strategy <- match.arg(start_strategy)
  base <- prep$lambda_frem
  base[!is.finite(base) | base <= 0] <- 1
  starts <- .fgee_fastk_start_grid(
    prep, qreml_lambda = qreml_lambda, strategy = start_strategy, bound = bound
  )

  score_u <- function(u, grad = FALSE) {
    lambda <- base * 10^as.numeric(u)
    fgee_fastk_score_grad(prep, lambda, need_gradient = grad)
  }
  start_score <- vapply(seq_len(nrow(starts)), function(i) {
    score_u(starts[i, ], FALSE)$score
  }, numeric(1))
  ord <- order(start_score)
  chosen <- ord[seq_len(min(as.integer(n_starts), length(ord)))]

  n_eval <- 0L
  opts <- vector("list", length(chosen))
  for (s in seq_along(chosen)) {
    last_u <- last <- NULL
    eval_both <- function(u) {
      uu <- as.numeric(u)
      if (!is.null(last_u) && identical(uu, last_u)) return(last)
      last_u <<- uu
      last <<- score_u(uu, TRUE)
      n_eval <<- n_eval + 1L
      last
    }
    fn <- function(u) {
      val <- eval_both(u)$score
      if (is.finite(val)) val else .Machine$double.xmax^0.25
    }
    gr <- function(u) {
      val <- eval_both(u)$gradient
      val[!is.finite(val)] <- 0
      val
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
  if (verbose) {
    message("gradient fastK score=", signif(best$value, 8),
            "; evaluations=", nrow(starts) + n_eval)
  }
  list(
    method = "fastk_grad",
    lambda = as.numeric(lambda),
    lambda.star = matrix(as.numeric(lambda), nrow = 1L),
    score = as.numeric(best$value),
    penalty_mat = penalty_from_setup(prep$ps, lambda),
    evaluations = nrow(starts) + n_eval,
    probe_evaluations = nrow(starts),
    gradient_evaluations = n_eval,
    convergence = best$convergence,
    message = best$message,
    boundary = any(abs(best$par) > bound - 0.05),
    log10_multiplier = as.numeric(best$par),
    optim = best,
    workspace = prep
  )
}
