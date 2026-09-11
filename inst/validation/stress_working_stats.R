suppressPackageStartupMessages({
  library(fastFGEE)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1L]))) else getwd()
source(file.path(script_dir, "helpers.R"))

cat("Running working-statistics stress checks\n")
seeds <- c(11L, 23L, 47L)
for (seed in seeds) {
  for (nl in c(1L, 2L, 7L)) {
    for (nf in c(1L, 3L, 11L)) {
      fx <- validation_fixture(N = 3L, n_long = nl, n_fun = nf, p = 5L,
                               rho_long = 0.7, rho_fun = -0.6, seed = seed + nl + nf)
      for (cl in c("independent", "ar1", "exchangeable")) {
        for (cf in c("independent", "ar1", "exchangeable")) {
          # Exchangeable rho must satisfy the dimension-dependent lower bound.
          if (cf == "exchangeable") fx$data[, rho_fn := 0.6]
          if (cl == "exchangeable") fx$data[, rho_long := 0.6]
          a <- fastFGEE:::fgee_build_working_stats(
            fx$data, fx$namesd, corr_fn = cf, corr_long = cl,
            retain = "scores", corr_solver = "exact"
          )
          b <- fastFGEE:::fgee_build_working_stats(
            fx$data[sample(.N)], fx$namesd, corr_fn = cf, corr_long = cl,
            retain = "scores", corr_solver = "exact"
          )
          stopifnot(isTRUE(all.equal(a$W_sum, b$W_sum, tolerance = 1e-9)))
          stopifnot(isTRUE(all.equal(a$D, b$D, tolerance = 1e-9)))
        }
      }
    }
  }
}

bad <- validation_fixture(N = 1L, n_long = 3L, n_fun = 4L, p = 3L)
err1 <- try(fastFGEE:::fgee_build_working_stats(
  bad$data[-1L], bad$namesd, corr_fn = "ar1", corr_long = "ar1"
), silent = TRUE)
err2 <- try(fastFGEE:::fgee_build_working_stats(
  rbind(bad$data, bad$data[1L]), bad$namesd,
  corr_fn = "ar1", corr_long = "ar1"
), silent = TRUE)
stopifnot(inherits(err1, "try-error"), inherits(err2, "try-error"))

cat("All working-statistics stress checks passed.\n")
