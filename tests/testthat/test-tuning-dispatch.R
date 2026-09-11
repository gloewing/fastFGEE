test_that("automatic selector dispatch is family specific", {
  expect_equal(
    fastFGEE:::.fgee_resolve_sp_method("auto", gaussian(), "identity"),
    "fastk_staged"
  )
  # Since 0.3.0 supported non-Gaussian families resolve to the fast
  # initialiser, which optimises the same exact fastK criterion from a much
  # cheaper start search. Validated at 300 replicates per configuration against
  # "fastk_grad" for binomial, Poisson and Gamma, under correct working
  # correlation and under three misspecification scenarios.
  expect_equal(
    fastFGEE:::.fgee_resolve_sp_method("auto", binomial(), "logit"),
    "fastk_grad_fast"
  )
  expect_equal(
    fastFGEE:::.fgee_resolve_sp_method("auto", poisson(), "log"),
    "fastk_grad_fast"
  )
  expect_equal(
    fastFGEE:::.fgee_resolve_sp_method("auto", Gamma(link = "log"), "log"),
    "fastk_grad_fast"
  )
})

test_that("qREML-fastK always returns an exact fastK criterion solution", {
  fx <- make_working_fixture(N = 18L, n_long = 2L, n_fun = 4L, p = 7L, seed = 131)
  fit <- make_fake_fit(p = 7L, q = 2L, family = binomial())
  fx$data[, Y := rbinom(.N, 1L, 0.4)]
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores"
  )

  ans <- fastFGEE:::fgee_tune_smoothing_optimized(
    working = st,
    dx = fx$data,
    namesd = fx$namesd,
    fit.initial = fit,
    cv.grid = list(c(0.1, 1, 10)),
    sp.method = "qreml_fastk",
    K = 3L,
    seed = 7L,
    fastk.start = "compact",
    keep.workspace = TRUE,
    verbose = FALSE
  )

  expect_equal(ans$method, "qreml_fastk")
  expect_true(is.list(ans$diagnostics$qreml))
  expect_true(is.list(ans$diagnostics$fastk))

  fprep <- ans$workspace$fastk
  qlambda <- ans$raw$qreml$lambda
  qscore <- fastFGEE:::fgee_fastk_score_grad(
    fprep, qlambda, need_gradient = FALSE
  )$score
  final_score <- fastFGEE:::fgee_fastk_score_grad(
    fprep, ans$lambda, need_gradient = FALSE
  )$score

  expect_lte(final_score, qscore + 1e-10)
  expect_equal(final_score, ans$score, tolerance = 1e-10)
})

test_that("compact tuning results do not retain workspaces by default", {
  fx <- make_working_fixture(N = 12L, n_long = 2L, n_fun = 4L, p = 7L, seed = 132)
  fit <- make_fake_fit(p = 7L, q = 2L, family = gaussian())
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores", exact_gaussian = TRUE,
    beta0 = fit$coefficients
  )
  prep <- fastFGEE:::fgee_fastk_prepare(
    st, fx$data, fx$namesd, fit, K = 3L, seed = 3L, exact = TRUE
  )
  raw <- fastFGEE:::fgee_tune_fastk_grad(
    prep, start_strategy = "compact", n_starts = 1L, maxit = 2L
  )
  compact <- fastFGEE:::.fgee_tune_result(raw, keep.workspace = FALSE)

  expect_null(compact$workspace)
  expect_null(compact$raw)
  expect_true(is.numeric(compact$lambda))
})

test_that("current public defaults use wild CIs and compact score retention", {
  expect_identical(eval(formals(fastFGEE::fgee)$joint.CI), "wild")
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      requested = "auto",
      joint.CI = "wild",
      var.type = "sandwich",
      sp.method = "fastk_staged"
    ),
    "scores"
  )
})

test_that("pure qREML uses aggregate retention while qREML-fastK needs scores", {
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      requested = "auto", joint.CI = FALSE, var.type = "sandwich",
      sp.method = "sandwich_qreml", gee.fit = TRUE
    ),
    "aggregate"
  )
  expect_equal(
    fastFGEE:::.fgee_resolve_working_retain(
      requested = "auto", joint.CI = FALSE, var.type = "sandwich",
      sp.method = "qreml_fastk", gee.fit = TRUE
    ),
    "scores"
  )
})
