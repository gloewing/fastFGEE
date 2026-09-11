validation_fixture <- function(N = 6L, n_long = 4L, n_fun = 8L, p = 7L,
                               rho_long = 0.35, rho_fun = 0.65,
                               family = c("gaussian", "binomial", "poisson", "gamma"),
                               seed = 481L) {
  family <- match.arg(family)
  set.seed(seed)
  rows <- as.integer(N * n_long * n_fun)
  dd <- data.table::CJ(
    cname_ = seq_len(N),
    time = seq_len(n_long),
    yindex.vec = seq_len(n_fun)
  )
  data.table::setorderv(dd, c("cname_", "time", "yindex.vec"))
  namesd <- paste0("X", seq_len(p))
  for (j in seq_len(p)) data.table::set(dd, j = namesd[j], value = stats::rnorm(rows))
  data.table::set(dd, j = namesd[1L], value = rep.int(1, rows))
  eta <- 0.15 + 0.1 * dd[[namesd[min(2L, p)]]] -
    0.08 * dd[[namesd[min(3L, p)]]]

  if (family == "gaussian") {
    dd[, Y := eta + stats::rnorm(.N)]
    dd[, v := exp(0.05 * time - 0.02 * yindex.vec)]
    dd[, muprime := 1]
    mu <- eta
  } else if (family == "binomial") {
    mu <- stats::plogis(eta)
    dd[, Y := stats::rbinom(.N, 1L, mu)]
    dd[, v := pmax(mu * (1 - mu), 1e-8)]
    dd[, muprime := v]
  } else if (family == "poisson") {
    mu <- exp(eta)
    dd[, Y := stats::rpois(.N, mu)]
    dd[, v := pmax(mu, 1e-8)]
    dd[, muprime := mu]
  } else {
    mu <- exp(eta)
    shape <- 3
    dd[, Y := stats::rgamma(.N, shape = shape, scale = mu / shape)]
    dd[, v := pmax(mu^2, 1e-8)]
    dd[, muprime := mu]
  }
  dd[, sqrtv := sqrt(v)]
  dd[, resid := (Y - mu) / sqrtv]
  dd[, rho_long := rho_long]
  dd[, rho_fn := rho_fun]
  list(data = dd, namesd = namesd)
}

validation_fake_fit <- function(p = 7L, nsdf = 1L, q = 2L,
                                family = stats::gaussian()) {
  stopifnot(p > nsdf + q)
  sizes <- rep((p - nsdf) %/% q, q)
  rem <- (p - nsdf) %% q
  if (rem > 0L) sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1L
  starts <- nsdf + cumsum(c(1L, head(sizes, -1L)))
  smooth <- lapply(seq_len(q), function(j) {
    k <- sizes[j]
    D <- diff(diag(k), differences = min(2L, k - 1L))
    list(
      first.para = starts[j],
      last.para = starts[j] + k - 1L,
      S = list(crossprod(D)),
      by = "NA"
    )
  })
  list(
    coefficients = seq_len(p) / 20,
    nsdf = nsdf,
    smooth = smooth,
    sp = seq_len(q),
    family = family
  )
}

validation_finite_diff <- function(fn, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (fn(xp) - fn(xm)) / (2 * eps)
  }, numeric(1))
}

validation_time <- function(fun, repeats = 3L) {
  if (!is.function(fun)) stop("fun must be a zero-argument function.")
  repeats <- as.integer(repeats)
  if (!is.finite(repeats) || repeats < 1L) stop("repeats must be positive.")
  elapsed <- numeric(repeats)
  value <- NULL
  for (i in seq_len(repeats)) {
    gc()
    tm <- system.time(value <- fun())
    elapsed[i] <- unname(tm[["elapsed"]])
  }
  list(value = value, median_elapsed = stats::median(elapsed), elapsed = elapsed)
}

validation_allocated_bytes <- function(fun) {
  if (!is.function(fun)) stop("fun must be a zero-argument function.")
  f <- tempfile("fgee-rprofmem-")
  on.exit(unlink(f), add = TRUE)
  utils::Rprofmem(f)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  value <- fun()
  utils::Rprofmem(NULL)
  lines <- readLines(f, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  list(value = value, bytes = sum(bytes, na.rm = TRUE), events = sum(is.finite(bytes)))
}
