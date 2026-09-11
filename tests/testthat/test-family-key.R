# Canonical family keys.
#
# The strings that matter are the ones a *fitted* mgcv extended family carries,
# because they contain the estimated parameter. An unfitted `mgcv::nb()` reports
# the clean "negative binomial" and its `getTheta(TRUE)` silently returns 1, so a
# fixture built from the constructor would pass a broken implementation while
# exercising neither the parsing nor a real nuisance value. Hence the literal
# strings below.

fam_obj <- function(family, link) {
  structure(list(family = family, link = link), class = "family")
}

test_that("fitted negative-binomial family strings resolve to one key", {
  for (s in c("Negative Binomial(2.338)", "Negative Binomial(2)",
              "Negative Binomial(13.227)", "negative binomial",
              "Negative Binomial", "negbin", "nbinom")) {
    expect_identical(fastFGEE:::.fgee_family_key(fam_obj(s, "log")),
                     "negbinomial", info = s)
  }
})

test_that("fitted beta-regression family strings resolve to one key", {
  for (s in c("Beta regression(13.227)", "Beta regression(4.399)",
              "Beta regression", "betar", "beta")) {
    expect_identical(fastFGEE:::.fgee_family_key(fam_obj(s, "logit")),
                     "beta", info = s)
  }
})

test_that("beta detection does not swallow the binomial families", {
  expect_identical(fastFGEE:::.fgee_family_key(stats::binomial()), "binomial")
  expect_identical(fastFGEE:::.fgee_family_key(fam_obj("quasibinomial", "logit")),
                   "quasibinomial")
  # A beta-binomial must not be claimed by the beta branch.
  expect_false(fastFGEE:::.fgee_family_key(fam_obj("Beta-Binomial", "logit")) ==
                 "beta")
})

test_that("the existing four families keep the keys they had", {
  expect_identical(fastFGEE:::.fgee_family_key(stats::gaussian()), "gaussian")
  expect_identical(fastFGEE:::.fgee_family_key(stats::poisson()), "poisson")
  expect_identical(fastFGEE:::.fgee_family_key(stats::Gamma(link = "log")),
                   "gamma")
  expect_identical(fastFGEE:::.fgee_family_key("Gamma"), "gamma")
  # quasi-variants stay distinct: delegating gamma would have collapsed this one
  expect_identical(fastFGEE:::.fgee_family_key(fam_obj("quasigamma", "log")),
                   "quasigamma")
  expect_identical(fastFGEE:::.fgee_family_key(fam_obj("quasipoisson", "log")),
                   "quasipoisson")
})
