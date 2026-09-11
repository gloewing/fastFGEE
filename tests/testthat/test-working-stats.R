test_that("one-pass W and D match legacy implementation", {
  fx <- make_working_fixture()
  combos <- expand.grid(
    corr_long = c("independent", "ar1", "exchangeable"),
    corr_fn = c("independent", "ar1", "exchangeable"),
    stringsAsFactors = FALSE
  )

  for (r in seq_len(nrow(combos))) {
    cl <- combos$corr_long[r]
    cf <- combos$corr_fn[r]
    oldW <- fastFGEE:::.getW(
      fx$data, fx$namesd, "cname_", corr_fn = cf, corr_long = cl,
      ensure_order = "setorderv"
    )
    oldD <- fastFGEE:::.getD(
      fx$data, fx$namesd, "cname_", corr_fn = cf, corr_long = cl,
      resid_col = "resid", ensure_order = "setorderv"
    )
    nw <- fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, corr_fn = cf, corr_long = cl,
      retain = "full", corr_solver = "exact"
    )
    expect_equal(lapply(nw$W, unname), lapply(oldW, unname), tolerance = 2e-8, info = paste(cl, cf))
    expect_equal(fastFGEE:::.fgee_working_dlist(nw), oldD,
                 tolerance = 2e-8, info = paste(cl, cf))
  }
})

test_that("row shuffling is corrected by canonical ordering", {
  fx <- make_working_fixture(seed = 17)
  a <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "scores", corr_solver = "exact"
  )
  set.seed(18)
  b <- fastFGEE:::fgee_build_working_stats(
    fx$data[sample(.N)], fx$namesd,
    corr_fn = "ar1", corr_long = "exchangeable",
    retain = "scores", corr_solver = "exact"
  )
  expect_equal(a$W_sum, b$W_sum, tolerance = 1e-10)
  expect_equal(a$D, b$D, tolerance = 1e-10)
})

test_that("exact Gaussian RHS identity is respected", {
  fx <- make_working_fixture(seed = 22)
  beta0 <- c(0.2, -0.4, 0.7)
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "full", corr_solver = "exact", exact_gaussian = TRUE,
    beta0 = beta0
  )
  for (i in seq_len(st$N)) {
    expect_equal(
      st$D_exact[, i],
      st$D[, i] + as.numeric(st$W[[i]] %*% beta0),
      tolerance = 1e-12
    )
  }
  expect_equal(unname(st$d_exact_sum), unname(rowSums(st$D_exact)), tolerance = 1e-12)
})

test_that("retention profiles preserve aggregates and reject unsafe choices", {
  fx <- make_working_fixture(seed = 33)
  a <- fastFGEE:::fgee_build_working_stats(fx$data, fx$namesd, retain = "aggregate")
  s <- fastFGEE:::fgee_build_working_stats(fx$data, fx$namesd, retain = "scores")
  f <- fastFGEE:::fgee_build_working_stats(fx$data, fx$namesd, retain = "full")
  expect_null(a$D); expect_null(a$W)
  expect_true(is.matrix(s$D)); expect_null(s$W)
  expect_true(is.matrix(f$D)); expect_true(is.list(f$W))
  expect_equal(a$W_sum, s$W_sum)
  expect_equal(s$W_sum, f$W_sum)
  expect_error(
    fastFGEE:::.fgee_resolve_working_retain(
      "aggregate", joint.CI = "wild", var.type = "sandwich"
    ),
    "insufficient"
  )
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      "auto", joint.CI = "wild", var.type = "sandwich"
    ),
    "scores"
  )
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      "auto", joint.CI = FALSE, var.type = "boot"
    ),
    "full"
  )
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      "auto", joint.CI = FALSE, var.type = "sandwich",
      sp.method = "sandwich_qreml"
    ),
    "aggregate"
  )
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      "auto", joint.CI = "wild", var.type = "sandwich", exact = TRUE
    ),
    "scores"
  )
})

test_that("incomplete and duplicate Kronecker grids fail loudly", {
  fx <- make_working_fixture(N = 1L)
  expect_error(
    fastFGEE:::fgee_build_working_stats(
      fx$data[-1L], fx$namesd, corr_fn = "ar1", corr_long = "ar1"
    ),
    "complete tensor-product grid"
  )
  dup <- rbind(fx$data, fx$data[1L])
  expect_error(
    fastFGEE:::fgee_build_working_stats(
      dup, fx$namesd, corr_fn = "ar1", corr_long = "ar1"
    ),
    "exactly one row"
  )
})

test_that("clusters may have unequal complete tensor-product grids", {
  a <- make_working_fixture(N = 1L, n_long = 3L, n_fun = 5L, p = 3L, seed = 44)
  b <- make_working_fixture(N = 1L, n_long = 5L, n_fun = 5L, p = 3L, seed = 45)
  b$data[, cname_ := 2L]
  dd <- data.table::rbindlist(list(a$data, b$data), use.names = TRUE)
  namesd <- a$namesd

  oldW <- fastFGEE:::.getW(
    dd, namesd, "cname_", corr_fn = "ar1", corr_long = "exchangeable",
    ensure_order = "setorderv"
  )
  oldD <- fastFGEE:::.getD(
    dd, namesd, "cname_", corr_fn = "ar1", corr_long = "exchangeable",
    resid_col = "resid", ensure_order = "setorderv"
  )
  nw <- fastFGEE:::fgee_build_working_stats(
    dd, namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "full", corr_solver = "exact"
  )

  expect_equal(nw$cluster_size, c(15L, 25L))
  expect_equal(lapply(nw$W, unname), lapply(oldW, unname), tolerance = 2e-8)
  expect_equal(fastFGEE:::.fgee_working_dlist(nw), oldD, tolerance = 2e-8)
})

test_that("scores retention avoids duplicated legacy cluster-list output", {
  fx <- make_working_fixture(seed = 52)
  s <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores"
  )
  f <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "full"
  )

  expect_null(fastFGEE:::.fgee_compat_scores(s))
  expect_null(fastFGEE:::.fgee_compat_bread(s))
  expect_true(is.list(fastFGEE:::.fgee_compat_scores(f)))
  expect_true(is.list(fastFGEE:::.fgee_compat_bread(f)))
})

test_that("one-direction batching and varying-rho fallbacks match legacy", {
  fx <- make_working_fixture(N = 2L, n_long = 5L, n_fun = 7L, p = 4L,
                             seed = 202)

  # Constant rho uses the batched one-direction path.
  for (spec in list(
    list(cf = "independent", cl = "ar1"),
    list(cf = "exchangeable", cl = "independent")
  )) {
    oldW <- fastFGEE:::.getW(
      fx$data, fx$namesd, "cname_", corr_fn = spec$cf,
      corr_long = spec$cl, ensure_order = "setorderv"
    )
    oldD <- fastFGEE:::.getD(
      fx$data, fx$namesd, "cname_", corr_fn = spec$cf,
      corr_long = spec$cl, resid_col = "resid",
      ensure_order = "setorderv"
    )
    nw <- fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, corr_fn = spec$cf, corr_long = spec$cl,
      retain = "full", corr_solver = "exact"
    )
    expect_equal(lapply(nw$W, unname), lapply(oldW, unname), tolerance = 2e-8)
    expect_equal(fastFGEE:::.fgee_working_dlist(nw), oldD,
                 tolerance = 2e-8)
  }

  # Longitudinal rho varying over the functional index must use the grouped
  # fallback and still reproduce the established implementation.
  fx$data[, rho_long := 0.15 + 0.6 *
            (yindex.vec - min(yindex.vec)) /
            max(1, max(yindex.vec) - min(yindex.vec))]
  oldW <- fastFGEE:::.getW(
    fx$data, fx$namesd, "cname_", corr_fn = "independent",
    corr_long = "ar1", ensure_order = "setorderv"
  )
  oldD <- fastFGEE:::.getD(
    fx$data, fx$namesd, "cname_", corr_fn = "independent",
    corr_long = "ar1", resid_col = "resid",
    ensure_order = "setorderv"
  )
  nw <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "independent", corr_long = "ar1",
    retain = "full", corr_solver = "exact"
  )
  expect_equal(lapply(nw$W, unname), lapply(oldW, unname), tolerance = 2e-8)
  expect_equal(fastFGEE:::.fgee_working_dlist(nw), oldD,
               tolerance = 2e-8)
})

test_that("working-statistics construction does not change data.table threads", {
  old <- data.table::getDTthreads()
  fx <- make_working_fixture(seed = 222)
  invisible(fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "scores", corr_solver = "exact"
  ))
  expect_identical(data.table::getDTthreads(), old)
})

test_that("gee.fit FALSE can use aggregate retention without wild inference", {
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      requested = "auto", joint.CI = FALSE, var.type = "sandwich",
      sp.method = "fastk_staged", gee.fit = FALSE
    ),
    "aggregate"
  )
})

test_that("full retention preserves legacy compatibility fields", {
  fx <- make_working_fixture(seed = 303)
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd,
    corr_fn = "ar1", corr_long = "exchangeable",
    retain = "full", corr_solver = "exact"
  )

  dl <- fastFGEE:::.fgee_compat_scores(st)
  wl <- fastFGEE:::.fgee_compat_bread(st)
  expect_length(dl, st$N)
  expect_length(wl, st$N)
  expect_equal(unname(do.call(cbind, dl)), unname(st$D), tolerance = 0)
  expect_equal(Reduce(`+`, wl), st$W_sum, tolerance = 1e-12)
})

test_that("unsorted data fail when canonical sorting is disabled", {
  fx <- make_working_fixture(seed = 304)
  bad <- data.table::copy(fx$data)
  bad <- bad[c(2L, 1L, seq.int(3L, nrow(bad)))]

  expect_error(
    fastFGEE:::fgee_build_working_stats(
      bad, fx$namesd,
      corr_fn = "ar1", corr_long = "exchangeable",
      retain = "scores", corr_solver = "exact",
      ensure_order = FALSE
    ),
    "canonical"
  )
})
