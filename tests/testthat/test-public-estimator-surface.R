test_that("public fgee exposes only the one-step estimator surface", {
  public_names <- names(formals(fastFGEE::fgee))
  expect_false(any(c(
    "exact", "gee.fit", "max.iter", "tune.method", "working.engine"
  ) %in% public_names))
  expect_true(all(c("sp.method", "fastk.memory", "corr.solver") %in%
                    public_names))

  sp_choices <- eval(formals(fastFGEE::fgee)$sp.method)
  expect_identical(
    sp_choices,
    c("auto", "fastk_staged", "fastk_grad", "fastk_grad_fast",
      "sandwich_qreml", "qreml_fastk")
  )
  expect_false("legacy" %in% sp_choices)

  expect_error(
    fastFGEE::fgee(exact = TRUE),
    "Only the validated one-step estimator"
  )
  expect_error(
    fastFGEE::fgee(max.iter = 10),
    "Only the validated one-step estimator"
  )
  expect_error(
    fastFGEE::fgee(working.engine = "legacy"),
    "Only the validated one-step estimator"
  )
  expect_error(
    fastFGEE::fgee(cv = "fullkfold"),
    "supports only cv"
  )
  expect_error(
    fastFGEE::fgee(unused_argument = 1),
    "does not silently ignore"
  )

  expect_true(is.function(fastFGEE:::.fgee_fit_internal))
  expect_true(is.function(fastFGEE:::.fgee_fit_internal_core))
  expect_true(is.function(fastFGEE:::.fgee_attach_nuisance_metadata))
  expect_false(".fgee_fit_internal" %in% getNamespaceExports("fastFGEE"))
})

test_that("auto dispatch remains staged Gaussian and fast gradient otherwise", {
  expect_identical(
    fastFGEE:::.fgee_resolve_sp_method("auto", gaussian(), "identity"),
    "fastk_staged"
  )
  expect_identical(
    fastFGEE:::.fgee_resolve_sp_method("auto", binomial(), "logit"),
    "fastk_grad_fast"
  )
  expect_identical(
    fastFGEE:::.fgee_resolve_sp_method("auto", mgcv::nb(), "log"),
    "fastk_grad_fast"
  )
  expect_identical(
    fastFGEE:::.fgee_resolve_sp_method("auto", mgcv::betar(), "logit"),
    "fastk_grad_fast"
  )
})

test_that("lowmem is explicitly deprecated before fitting begins", {
  result <- NULL
  expect_warning(
    result <- try(fastFGEE::fgee(fastk.memory = "lowmem"), silent = TRUE),
    "lowmem.*deprecated"
  )
  expect_s3_class(result, "try-error")
})
