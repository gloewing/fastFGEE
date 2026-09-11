test_that("Gamma profile and moment estimators recover positive dispersion", {
  set.seed(3101)
  n <- 12000L
  mu <- exp(seq(log(0.5), log(3), length.out = n))
  phi <- 0.45
  y <- stats::rgamma(n, shape = 1 / phi, scale = phi * mu)
  pr <- fastFGEE:::fgee_estimate_nuisance(y, mu, stats::Gamma(link = "log"), method = "profile")
  mo <- fastFGEE:::fgee_estimate_nuisance(y, mu, stats::Gamma(link = "log"), method = "moment")
  expect_equal(pr$parameter, "dispersion")
  expect_lt(abs(pr$value - phi), 0.06)
  expect_lt(abs(mo$value - phi), 0.07)
})

test_that("beta precision profile and stable moment estimator behave sensibly", {
  set.seed(3102)
  n <- 15000L
  mu <- stats::plogis(seq(-2, 2, length.out = n))
  precision <- 24
  y <- stats::rbeta(n, mu * precision, (1 - mu) * precision)
  fam <- structure(list(family = "Beta regression", link = "logit"), class = "family")
  pr <- fastFGEE:::fgee_estimate_nuisance(y, mu, fam, method = "profile")
  mo <- fastFGEE:::fgee_estimate_nuisance(y, mu, fam, method = "moment")
  expect_equal(pr$parameter, "precision")
  expect_lt(abs(pr$value - precision) / precision, 0.12)
  expect_lt(abs(mo$value - precision) / precision, 0.18)
})

test_that("negative-binomial theta profile and moment estimator are valid", {
  set.seed(3103)
  n <- 15000L
  mu <- exp(seq(log(0.3), log(5), length.out = n))
  theta <- 5
  y <- stats::rnbinom(n, size = theta, mu = mu)
  fam <- structure(list(family = "Negative Binomial", link = "log"), class = "family")
  pr <- fastFGEE:::fgee_estimate_nuisance(y, mu, fam, method = "profile")
  mo <- fastFGEE:::fgee_estimate_nuisance(y, mu, fam, method = "moment")
  expect_equal(pr$parameter, "theta")
  expect_lt(abs(pr$value - theta) / theta, 0.18)
  expect_lt(abs(mo$value - theta) / theta, 0.25)
})

test_that("negative-binomial underdispersion is reported as Poisson boundary", {
  mu <- rep(2, 500)
  y <- rep(c(1, 2, 3, 2), length.out = 500)
  fam <- structure(list(family = "Negative Binomial", link = "log"), class = "family")
  mo <- fastFGEE:::fgee_estimate_nuisance(y, mu, fam, method = "moment")
  expect_true(mo$boundary)
  expect_gt(mo$value, 1e6)
})

test_that("explicit fixed nuisance is retained and variance columns are consistent", {
  dd <- data.table::data.table(Y = c(0.6, 1.1, 2.2, 0.9), x = c(-1, 0, 1, 2))
  ans <- fastFGEE:::fgee_update_working_cols_dt(
    dx = dd, namesd = "x", beta = 0.2,
    family = stats::Gamma(link = "log"), link = "log",
    update_nuisance = "fixed", dispersion = 0.4, copy = TRUE
  )
  n <- attr(ans, "nuisance")
  expect_equal(n$dispersion, 0.4)
  expect_equal(ans$v, 0.4 * ans$p^2, tolerance = 1e-13)
  expect_equal(ans$resid, (ans$Y - ans$p) / sqrt(ans$v), tolerance = 1e-13)
})

test_that("missing Gamma dispersion is estimated with an explicit warning", {
  set.seed(3104)
  dd <- data.table::data.table(Y = stats::rgamma(300, shape = 2, scale = 0.5), x = stats::rnorm(300))
  expect_warning(
    ans <- fastFGEE:::fgee_update_working_cols_dt(
      dx = dd, namesd = "x", beta = 0,
      family = stats::Gamma(link = "log"), link = "log",
      update_nuisance = "fixed", copy = TRUE
    ),
    "rather than silently using 1"
  )
  n <- attr(ans, "nuisance")
  expect_s3_class(n, "fgee_nuisance")
  expect_true(is.finite(n$dispersion) && n$dispersion > 0)
})

test_that("beta boundary handling is explicit and recorded", {
  fam <- structure(list(family = "Beta regression", link = "logit"), class = "family")
  expect_warning(
    ans <- fastFGEE:::fgee_estimate_nuisance(
      y = c(0, 0.2, 0.7, 1), mu = c(0.1, 0.3, 0.6, 0.9),
      family = fam, method = "moment"
    ),
    "clipped"
  )
  expect_equal(ans$clipped_y, 2L)
})

test_that("working-statistic objects preserve typed nuisance metadata", {
  dd <- data.table::data.table(
    cname_ = rep(1:2, each = 4),
    time = rep(rep(1:2, each = 2), 2),
    yindex.vec = rep(1:2, 4),
    Y = c(0.8, 1.2, 1.1, 0.7, 1.4, 0.9, 1.3, 0.6),
    X1 = 1,
    p = 1,
    muprime = 1,
    v = 0.5,
    sqrtv = sqrt(0.5),
    resid = 0,
    rho_fn = 0,
    rho_long = 0
  )
  data.table::setattr(dd, "nuisance", structure(
    list(family = "gamma", parameter = "dispersion", value = 0.5,
         dispersion = 0.5, method = "fixed", converged = TRUE,
         boundary = FALSE),
    class = c("fgee_nuisance", "list")
  ))
  st <- fastFGEE:::fgee_build_working_stats(
    dx = dd, namesd = "X1", cname_ = "cname_",
    corr_fn = "independent", corr_long = "independent",
    index_fn = "yindex.vec", index_long = "time",
    retain = "scores", copy_dt = TRUE
  )
  expect_equal(st$nuisance$dispersion, 0.5)
})
