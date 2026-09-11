make_working_fixture <- function(N = 3L, n_long = 4L, n_fun = 6L, p = 3L,
                                 rho_long = 0.35, rho_fun = 0.65,
                                 seed = 481) {
  set.seed(seed)
  rows <- N * n_long * n_fun
  dd <- data.table::CJ(
    cname_ = seq_len(N),
    time = seq_len(n_long),
    yindex.vec = seq_len(n_fun)
  )
  data.table::setorderv(dd, c("cname_", "time", "yindex.vec"))
  for (j in seq_len(p)) data.table::set(dd, j = paste0("X", j), value = stats::rnorm(rows))
  data.table::set(dd, j = "X1", value = rep.int(1, rows))
  dd[, v := exp(0.08 * time - 0.03 * yindex.vec)]
  dd[, sqrtv := sqrt(v)]
  dd[, muprime := 0.7 + 0.1 * abs(sin(yindex.vec))]
  dd[, resid := stats::rnorm(.N) / sqrtv]
  dd[, Y := stats::rnorm(.N)]
  dd[, rho_long := rho_long]
  dd[, rho_fn := rho_fun]
  list(data = dd, namesd = paste0("X", seq_len(p)))
}

make_fake_fit <- function(p = 7L, nsdf = 1L, q = 2L,
                          family = stats::gaussian()) {
  stopifnot(p > nsdf + q)
  sizes <- rep((p - nsdf) %/% q, q)
  sizes[seq_len((p - nsdf) %% q)] <- sizes[seq_len((p - nsdf) %% q)] + 1L
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

finite_diff <- function(fn, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (fn(xp) - fn(xm)) / (2 * eps)
  }, numeric(1))
}
