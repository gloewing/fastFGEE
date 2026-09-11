suppressPackageStartupMessages({
  library(fastFGEE)
  library(refund)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else getwd()
source(file.path(script_dir, "helpers.R"))

out_dir <- Sys.getenv("FGEE_VALIDATION_OUT", file.path(script_dir, "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
repeats <- as.integer(Sys.getenv("FGEE_BENCH_REPEATS", "2"))

data("d", package = "fastFGEE")
fit0 <- refund::pffr(
  Y ~ X1 + X2, data = d, family = binomial(link = "logit"),
  algorithm = "bam", discrete = TRUE,
  bs.yindex = list(bs = "ps", k = 7, m = c(2, 1)),
  bs.int = list(bs = "ps", k = 7, m = c(2, 1))
)
grids <- list(c(0.1, 1, 10), c(0.1, 1, 10), c(0.5, 1, 2))

run_legacy <- function() fastFGEE:::.fgee_fit_internal(
  Y ~ X1 + X2, data = d, cluster = "ID",
  family = binomial(link = "logit"), time = "time",
  corr_long = "exchangeable", corr_fn = "independent",
  pffr.mod = fit0, cv.grid = grids, joint.CI = FALSE,
  working.engine = "legacy", sp.method = "legacy", rho.smooth = FALSE
)
run_opt <- function() fastFGEE::fgee(
  Y ~ X1 + X2, data = d, cluster = "ID",
  family = binomial(link = "logit"), time = "time",
  corr_long = "exchangeable", corr_fn = "independent",
  pffr.mod = fit0, cv.grid = grids, joint.CI = FALSE,
  sp.method = "fastk_staged",
  working.retain = "scores", corr.solver = "supergauss",
  rho.smooth = FALSE, verbose.tuning = FALSE
)

legacy <- validation_time(run_legacy, repeats)
opt <- validation_time(run_opt, repeats)
stopifnot(max(abs(legacy$value$beta - opt$value$beta)) < 2e-7)
stopifnot(max(abs(legacy$value$vb - opt$value$vb)) < 5e-7)

ans <- data.frame(
  legacy_sec = legacy$median_elapsed,
  optimized_sec = opt$median_elapsed,
  speedup = legacy$median_elapsed / opt$median_elapsed,
  legacy_object_bytes = as.numeric(object.size(legacy$value)),
  optimized_object_bytes = as.numeric(object.size(opt$value)),
  beta_max_abs_error = max(abs(legacy$value$beta - opt$value$beta)),
  covariance_max_abs_error = max(abs(legacy$value$vb - opt$value$vb))
)
utils::write.csv(ans, file.path(out_dir, "fgee_engine_benchmark.csv"), row.names = FALSE)
print(ans, row.names = FALSE)
