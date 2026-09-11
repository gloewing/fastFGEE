suppressPackageStartupMessages({
  library(fastFGEE)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else getwd()
source(file.path(script_dir, "helpers.R"))

out_dir <- Sys.getenv("FGEE_VALIDATION_OUT", file.path(script_dir, "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
repeats <- as.integer(Sys.getenv("FGEE_BENCH_REPEATS", "3"))

spec <- data.frame(
  case = c("small", "many_visits", "many_clusters"),
  N = c(25L, 25L, 100L),
  n_long = c(5L, 100L, 5L),
  n_fun = c(50L, 50L, 50L),
  p = c(16L, 16L, 16L)
)
if (identical(Sys.getenv("FGEE_BENCH_LARGE", "0"), "1")) {
  spec <- rbind(spec, data.frame(
    case = "large", N = 100L, n_long = 100L, n_fun = 50L, p = 16L
  ))
}
custom_N <- Sys.getenv("FGEE_BENCH_N", "")
if (nzchar(custom_N)) {
  spec <- rbind(spec, data.frame(
    case = "custom",
    N = as.integer(custom_N),
    n_long = as.integer(Sys.getenv("FGEE_BENCH_NI", "100")),
    n_fun = as.integer(Sys.getenv("FGEE_BENCH_L", "50")),
    p = as.integer(Sys.getenv("FGEE_BENCH_P", "16"))
  ))
}

rows <- list()
for (z in seq_len(nrow(spec))) {
  s <- spec[z, ]
  fx <- validation_fixture(
    N = s$N, n_long = s$n_long, n_fun = s$n_fun, p = s$p,
    rho_long = 0.35, rho_fun = 0.75, seed = 1000L + z
  )

  old <- validation_time(function() {
    W <- fastFGEE:::.getW(
      fx$data, fx$namesd, "cname_", corr_fn = "ar1",
      corr_long = "exchangeable", ensure_order = "setorderv"
    )
    D <- fastFGEE:::.getD(
      fx$data, fx$namesd, "cname_", corr_fn = "ar1",
      corr_long = "exchangeable", resid_col = "resid",
      ensure_order = "setorderv"
    )
    list(W = W, D = D)
  }, repeats)

  new_exact <- validation_time(function() {
    fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
      retain = "full", corr_solver = "exact"
    )
  }, repeats)

  new_sg <- validation_time(function() {
    fastFGEE:::fgee_build_working_stats(
      fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
      retain = "full", corr_solver = "supergauss"
    )
  }, repeats)

  werr <- max(vapply(seq_along(old$value$W), function(i) {
    max(abs(old$value$W[[i]] - new_exact$value$W[[i]]))
  }, numeric(1)))
  derr <- max(vapply(seq_along(old$value$D), function(i) {
    max(abs(old$value$D[[i]] - new_exact$value$D[, i]))
  }, numeric(1)))
  sgerr <- max(abs(new_sg$value$W_sum - new_exact$value$W_sum),
               abs(new_sg$value$D - new_exact$value$D))
  if (max(werr, derr, sgerr) > 2e-7) stop("working-statistics mismatch")

  compact <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "scores", corr_solver = "exact"
  )
  aggregate <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, corr_fn = "ar1", corr_long = "exchangeable",
    retain = "aggregate", corr_solver = "exact"
  )

  rows[[length(rows) + 1L]] <- data.frame(
    case = s$case, N = s$N, n_long = s$n_long, n_fun = s$n_fun,
    p = s$p, M = s$N * s$n_long * s$n_fun,
    legacy_W_plus_D_sec = old$median_elapsed,
    one_pass_exact_sec = new_exact$median_elapsed,
    one_pass_supergauss_sec = new_sg$median_elapsed,
    legacy_over_exact_speedup = old$median_elapsed / new_exact$median_elapsed,
    full_bytes = as.numeric(object.size(new_exact$value)),
    scores_bytes = as.numeric(object.size(compact)),
    aggregate_bytes = as.numeric(object.size(aggregate)),
    max_abs_W_error = werr,
    max_abs_D_error = derr,
    max_abs_exact_supergauss_error = sgerr
  )
}

ans <- do.call(rbind, rows)
utils::write.csv(ans, file.path(out_dir, "working_stats_benchmark.csv"), row.names = FALSE)
print(ans, row.names = FALSE)
cat("Wrote", file.path(out_dir, "working_stats_benchmark.csv"), "\n")
