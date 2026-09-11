suppressPackageStartupMessages({
  library(fastFGEE)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else getwd()
source(file.path(script_dir, "helpers.R"))

out_dir <- Sys.getenv("FGEE_VALIDATION_OUT", file.path(script_dir, "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sizes <- as.integer(strsplit(Sys.getenv("FGEE_CORR_SIZES", "50,250,1000"), ",", fixed = TRUE)[[1L]])
repeats <- as.integer(Sys.getenv("FGEE_BENCH_REPEATS", "3"))
rhs_cols <- as.integer(Sys.getenv("FGEE_CORR_RHS", "24"))

rows <- list()
set.seed(901)
for (corr in c("ar1", "exchangeable")) {
  for (n in sizes) {
    rho <- if (corr == "ar1") 0.75 else 0.35
    B <- matrix(rnorm(n * rhs_cols), n, rhs_cols)

    exact <- validation_time(function() {
      fastFGEE:::.fgee_apply_corr_inverse(
        B, corr = corr, rho = rho, solver = "exact"
      )
    }, repeats)
    sg <- validation_time(function() {
      fastFGEE:::.fgee_apply_corr_inverse(
        B, corr = corr, rho = rho, solver = "supergauss"
      )
    }, repeats)

    err <- max(abs(exact$value - sg$value))
    if (!is.finite(err) || err > 2e-7) {
      stop(corr, " n=", n, ": exact and SuperGauss differ; max error=", err)
    }

    rows[[length(rows) + 1L]] <- data.frame(
      correlation = corr,
      n = n,
      rhs = rhs_cols,
      exact_sec = exact$median_elapsed,
      supergauss_sec = sg$median_elapsed,
      supergauss_over_exact = sg$median_elapsed / exact$median_elapsed,
      max_abs_error = err
    )
  }
}

ans <- do.call(rbind, rows)
utils::write.csv(ans, file.path(out_dir, "correlation_operator_benchmark.csv"), row.names = FALSE)
print(ans, row.names = FALSE)
cat("Wrote", file.path(out_dir, "correlation_operator_benchmark.csv"), "\n")
