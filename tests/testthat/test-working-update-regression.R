test_that("repeated working-column updates refresh muprime", {

  cases <- list(
    list(
      family = binomial(link = "logit"),
      y = c(0, 1, 0, 1, 1, 0)
    ),
    list(
      family = poisson(link = "log"),
      y = c(0, 1, 2, 1, 3, 2)
    ),
    list(
      family = Gamma(link = "log"),
      y = c(0.5, 1.0, 1.8, 0.9, 2.5, 1.4)
    )
  )

  for (z in cases) {

    d <- data.table::data.table(
      Y = z$y,
      x1 = c(-1.2, -0.4, 0.1, 0.7, 1.1, 1.8),
      x2 = c(0.3, 1.0, -0.5, 0.2, -1.1, 0.8)
    )

    a <- fastFGEE:::fgee_update_working_cols_dt(
      dx = d,
      namesd = c("x1", "x2"),
      beta = c(0.15, -0.20),
      family = z$family,
      link = z$family$link,
      exact = FALSE,
      update_nuisance = "fixed",
      copy = TRUE
    )

    mp_initial <- a$muprime

    b <- fastFGEE:::fgee_update_working_cols_dt(
      dx = a,
      namesd = c("x1", "x2"),
      beta = c(0.80, -0.55),
      family = z$family,
      link = z$family$link,
      exact = FALSE,
      update_nuisance = "fixed",
      copy = FALSE
    )

    expected <- fastFGEE:::.fgee_muprime_from_mu(
      b$p,
      link = z$family$link,
      clamp_eps = 1e-8
    )

    expect_equal(
      b$muprime,
      expected,
      tolerance = 1e-13,
      info = paste(z$family$family, z$family$link)
    )

    expect_gt(
      max(abs(b$muprime - mp_initial)),
      1e-6
    )
  }
})


test_that("repeated binomial update has muprime equal to variance function", {

  d <- data.table::data.table(
    Y = c(0, 1, 0, 1, 1, 0),
    x1 = seq(-1, 1, length.out = 6),
    x2 = c(1, -1, 0.5, -0.5, 0.2, -0.2)
  )

  a <- fastFGEE:::fgee_update_working_cols_dt(
    dx = d,
    namesd = c("x1", "x2"),
    beta = c(0.1, 0.2),
    family = binomial(link = "logit"),
    copy = TRUE
  )

  b <- fastFGEE:::fgee_update_working_cols_dt(
    dx = a,
    namesd = c("x1", "x2"),
    beta = c(0.9, -0.6),
    family = binomial(link = "logit"),
    copy = FALSE
  )

  expect_equal(
    b$muprime,
    b$p * (1 - b$p),
    tolerance = 1e-14
  )

  expect_equal(
    b$muprime,
    b$v,
    tolerance = 1e-14
  )
})
