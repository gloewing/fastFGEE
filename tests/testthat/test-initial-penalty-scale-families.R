# gee.fit = FALSE for dispersion families: Gamma and quasi-Poisson --------------
#
# Companion to test-initial-fit-sandwich.R, which covers Gaussian and Poisson.
# These two families discriminate the correct factor G * c_work from the two
# plausible wrong ones:
#
#   Gamma          c_work is the dispersion actually multiplied into the working
#                  variance (v = phi * mu^2), so the factor is G * phi, not G.
#   quasi-Poisson  the current fixed-nuisance path multiplies no extra dispersion
#                  into the working variance (v = mu), so the factor is G, not
#                  G * fit$sig2 -- the fixture deliberately has fitted sig2 far
#                  from one so that a G * fit$sig2 implementation would fail.
#
# The optimized and legacy engines may legitimately use DIFFERENT working
# dispersions for Gamma (the legacy path estimates one from the fitted mean
# rather than taking fit$sig2).  Each engine is therefore checked against the
# dispersion it actually used, established from its own recorded bread, and no
# equality is demanded between the two engines' covariances.
#
# Every expected value is assembled without calling the scaling, penalty or
# covariance code under test; those are called only to obtain the value tested.

ips_sim <- function(response, G = 25L, n_i = 3L, L = 12L, seed = 4L) {
  set.seed(seed)
  s <- seq(0, 1, length.out = L)
  ID <- rep(seq_len(G), each = n_i)
  n <- length(ID)
  X1 <- stats::rnorm(n)
  u <- matrix(stats::rnorm(2L * G), G, 2L)
  eta <- outer(rep(1, n), 0.4 + 0.4 * sin(2 * pi * s)) +
    X1 %o% (0.4 * cos(2 * pi * s)) +
    0.4 * u[ID, 1L] %o% rep(1, L) + 0.4 * u[ID, 2L] %o% s
  Y <- switch(response,
    gamma = matrix(stats::rgamma(n * L, shape = 2, rate = 2 / exp(eta)), n, L),
    poisson = matrix(stats::rpois(n * L, exp(eta)), n, L)
  )
  d <- data.frame(ID = ID, time = rep(seq_len(n_i), G), X1 = X1)
  d$Y <- I(Y)
  list(data = d, s = s, G = G, n_i = n_i, L = L)
}

ips_fit_pffr <- function(sim, family) {
  refund::pffr(
    Y ~ X1, yind = sim$s, data = sim$data, family = family,
    algorithm = "bam", discrete = TRUE, sandwich = "none",
    bs.yindex = list(bs = "ps", k = 6, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 6, m = c(2, 1))
  )
}

ips_run <- function(sim, fit, family, engine) {
  suppressMessages(suppressWarnings(fastFGEE:::.fgee_fit_internal(
    Y ~ X1, data = sim$data, cluster = "ID", family = family, time = "time",
    pffr.mod = fit, gee.fit = FALSE, var.type = "sandwich",
    corr_fn = "independent", corr_long = "independent",
    bs = "ps", knots = 5, joint.CI = FALSE, working.engine = engine,
    sp.method = if (identical(engine, "legacy")) "legacy" else "auto"
  )))
}

# Independent total-scale penalty from the fitted smooth components.  smooth$S is
# already in the fitted parameterisation, so it is NOT multiplied by S.scale; and
# S is never recovered from Vp, which would make the reference a function of the
# quantity under test.
ips_S_mgcv <- function(fit) {
  p <- length(fit$coefficients)
  S <- matrix(0, p, p)
  spv <- fit$full.sp
  if (is.null(spv) || !length(spv)) spv <- fit$sp
  spv <- as.numeric(spv)
  k <- 0L
  for (i in seq_along(fit$smooth)) {
    sm <- fit$smooth[[i]]
    idx <- sm$first.para:sm$last.para
    for (j in seq_along(sm$S)) {
      k <- k + 1L
      Sj <- as.matrix(sm$S[[j]])
      if (nrow(Sj) != length(idx)) {
        stop("penalty component ", k, ": S dim ", nrow(Sj),
             " != coefficient block ", length(idx))
      }
      S[idx, idx] <- S[idx, idx] + spv[k] * Sj
    }
  }
  if (k != length(spv)) {
    stop("assembled ", k, " penalty components but sp has ", length(spv))
  }
  list(S = S, n_comp = k, sp = spv)
}

# The engine's recorded bread and score.  The optimized path fills working0; the
# legacy path returns the per-cluster lists wi0 / di0 instead.
ips_recorded <- function(ans, engine) {
  if (identical(engine, "legacy")) {
    list(W_sum = Reduce(`+`, ans$wi0), d_sum = as.numeric(Reduce(`+`, ans$di0)),
         G = length(ans$wi0))
  } else {
    list(W_sum = as.matrix(ans$working0$W_sum),
         d_sum = as.numeric(ans$working0$d_sum), G = ans$working0$N)
  }
}

ips_fro <- function(x) sqrt(sum(as.matrix(x)^2))

ips_snapshot <- function(fit) {
  list(coef = fit$coefficients, sp = fit$sp, sig2 = fit$sig2, Vp = fit$Vp,
       edf = fit$edf, fitted = fit$fitted.values)
}

# One family case, run against both engines.
ips_check_family <- function(family, response, v_base_fun, expect_unit_phi,
                             stationarity_tol) {
  sim <- ips_sim(response)
  fit <- ips_fit_pffr(sim, family)
  before <- ips_snapshot(fit)

  # the fixture must be discriminating: sig2 clearly away from 1
  expect_gt(abs(as.numeric(fit$sig2)[1L] - 1), 0.1)

  ref_pen <- ips_S_mgcv(fit)
  expect_gt(ref_pen$n_comp, 1L)

  # independent state from the pffr fit, with the row key verified rather than
  # assumed: lpmatrix rows are curve-major over the response matrix
  X <- suppressWarnings(stats::predict(fit, type = "lpmatrix"))
  yv <- as.numeric(fit$y)
  mu <- as.numeric(fit$fitted.values)
  expect_equal(max(abs(as.numeric(t(as.matrix(sim$data$Y))) - yv)), 0)
  cl <- rep(sim$data$ID, each = sim$L)
  expect_identical(length(cl), nrow(X))
  expect_identical(colnames(X), names(fit$coefficients))

  mup <- mu                       # log link for both families here
  v_base <- v_base_fun(mu)
  W_unit <- crossprod(X * (mup / sqrt(v_base)))   # correct weighted Gram, phi = 1

  for (engine in c("optimized", "legacy")) {
    ans <- ips_run(sim, fit, family, engine)
    info <- paste(family$family, engine)
    rec <- ips_recorded(ans, engine)

    # 1. initial coefficients preserved, order checked explicitly
    expect_equal(unname(as.numeric(ans$beta)), unname(as.numeric(fit$coefficients)),
                 tolerance = 1e-10, info = info)
    expect_identical(
      intersect(names(as.data.frame(ans$data)), names(fit$coefficients)),
      names(fit$coefficients)
    )

    # 2. G counts independent subjects, not curves or scalar observations
    G <- rec$G
    expect_identical(as.integer(G), sim$G, info = info)
    expect_false(identical(as.integer(G), nrow(sim$data)))
    expect_false(identical(as.integer(G), nrow(X)))

    # 3. the working variance this engine actually used.  phi is recovered from
    #    the engine's own recorded bread and then confirmed by rebuilding that
    #    bread and score from scratch.
    phi <- sum(diag(W_unit)) / sum(diag(rec$W_sum))
    expect_true(is.finite(phi) && phi > 0, info = info)
    Z <- X * (mup / sqrt(phi * v_base))
    e <- (yv - mu) / sqrt(phi * v_base)
    expect_lt(ips_fro(crossprod(Z) - rec$W_sum) / ips_fro(rec$W_sum), 1e-10)
    u <- rowsum(Z * e, cl)
    expect_lt(ips_fro(colSums(u) - rec$d_sum) / ips_fro(rec$d_sum), 1e-10)

    if (expect_unit_phi) {
      # quasi path: no extra dispersion in the working variance
      expect_equal(phi, 1, tolerance = 1e-8, info = info)
    } else {
      expect_gt(abs(phi - 1), 0.1)      # Gamma: genuinely non-unit
    }

    # when the engine exposes the working columns, check v directly and against
    # its own nuisance record
    dd <- as.data.frame(ans$data)
    if (all(c("v", "p") %in% names(dd))) {
      expect_equal(as.numeric(dd$v), phi * v_base_fun(as.numeric(dd$p)),
                   tolerance = 1e-8, info = info)
      disp <- suppressWarnings(as.numeric(ans$working0$nuisance$dispersion)[1L])
      if (length(disp) && is.finite(disp)) {
        expect_equal(disp, phi, tolerance = 1e-8, info = info)
      }
    }

    c_work <- phi

    # 4/5. G * c_work * P_returned reproduces the independent total-scale penalty,
    #      and the returned lambda is the correctly mapped sp / (G * c_work)
    P <- as.matrix(ans$pen.mat)
    expect_lt(ips_fro(G * c_work * P - ref_pen$S) / ips_fro(ref_pen$S), 1e-10)
    expect_equal(as.numeric(ans$lambda) * G * c_work, ref_pen$sp,
                 tolerance = 1e-8, info = info)

    # 6. stationarity of the penalised working score at the initial fit
    db <- rec$d_sum / G
    expect_lt(
      ips_fro(db - as.numeric(P %*% as.numeric(fit$coefficients))) / ips_fro(db),
      stationarity_tol
    )

    # 7. covariance against the independent centred cluster sandwich
    H_sum <- crossprod(Z) + ref_pen$S / c_work
    B_sum <- solve(H_sum)
    uc <- sweep(u, 2L, colMeans(u))
    V_ref <- G / (G - 1) * B_sum %*% crossprod(uc) %*% t(B_sum)
    V <- unname(as.matrix(ans$vb))
    expect_lt(ips_fro(V - V_ref) / ips_fro(V_ref), 1e-6)
    expect_equal(unname(sqrt(diag(V))), unname(sqrt(diag(V_ref))),
                 tolerance = 1e-6, info = info)
    expect_lt(max(abs(V - t(V))) / max(abs(V)), 1e-10)
    expect_true(all(is.finite(V)), info = info)
  }

  # the shared initial fit was not modified in place by either engine
  expect_equal(ips_snapshot(fit), before, tolerance = 0)
}

test_that("gee.fit = FALSE penalty scale and sandwich are correct for Gamma", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")
  # Gamma has a non-canonical log link here, so stationarity of the penalised
  # working score holds to the IRLS convergence level rather than exactly.
  ips_check_family(stats::Gamma(link = "log"), "gamma",
                   v_base_fun = function(mu) mu^2,
                   expect_unit_phi = FALSE, stationarity_tol = 1e-5)
})

test_that("gee.fit = FALSE penalty scale and sandwich are correct for quasi-Poisson", {
  skip_on_cran()
  skip_if_not_installed("refund")
  skip_if_not_installed("mgcv")
  # Canonical log link with no extra working dispersion: stationarity is exact.
  ips_check_family(stats::quasipoisson(link = "log"), "poisson",
                   v_base_fun = function(mu) mu,
                   expect_unit_phi = TRUE, stationarity_tol = 1e-8)
})
