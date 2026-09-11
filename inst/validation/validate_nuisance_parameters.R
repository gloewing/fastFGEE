library(fastFGEE)

set.seed(93001)
results <- list()

n <- 20000L
mu_g <- exp(seq(log(0.4), log(4), length.out = n))
phi <- 0.35
y_g <- rgamma(n, shape = 1 / phi, scale = phi * mu_g)
results$gamma_profile <- fastFGEE:::fgee_estimate_nuisance(y_g, mu_g, Gamma(link="log"), "profile")
results$gamma_moment <- fastFGEE:::fgee_estimate_nuisance(y_g, mu_g, Gamma(link="log"), "moment")

mu_b <- plogis(seq(-2.5, 2.5, length.out = n))
prec <- 30
y_b <- rbeta(n, mu_b * prec, (1 - mu_b) * prec)
bfam <- structure(list(family="Beta regression", link="logit"), class="family")
results$beta_profile <- fastFGEE:::fgee_estimate_nuisance(y_b, mu_b, bfam, "profile")
results$beta_moment <- fastFGEE:::fgee_estimate_nuisance(y_b, mu_b, bfam, "moment")

mu_n <- exp(seq(log(0.2), log(6), length.out = n))
th <- 4
y_n <- rnbinom(n, size = th, mu = mu_n)
nfam <- structure(list(family="Negative Binomial", link="log"), class="family")
results$nb_profile <- fastFGEE:::fgee_estimate_nuisance(y_n, mu_n, nfam, "profile")
results$nb_moment <- fastFGEE:::fgee_estimate_nuisance(y_n, mu_n, nfam, "moment")

print(lapply(results, function(x) x[c("family","parameter","value","method","boundary","n_used")]))
stopifnot(abs(results$gamma_profile$value - phi) < 0.04)
stopifnot(abs(results$beta_profile$value - prec) / prec < 0.08)
stopifnot(abs(results$nb_profile$value - th) / th < 0.12)
cat("All nuisance-parameter validation checks passed.\n")
