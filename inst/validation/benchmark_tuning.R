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

family_object <- function(x) switch(
  x,
  gaussian = gaussian(),
  binomial = binomial(),
  poisson = poisson(),
  gamma = Gamma(link = "log")
)

families <- strsplit(
  Sys.getenv("FGEE_BENCH_FAMILIES", "gaussian,binomial,poisson,gamma"),
  ",", fixed = TRUE
)[[1L]]
spec <- expand.grid(
  family = families,
  case = c("small", "many_visits"),
  stringsAsFactors = FALSE
)
if (nzchar(Sys.getenv("FGEE_BENCH_N", ""))) {
  spec <- rbind(spec, expand.grid(
    family = families, case = "custom", stringsAsFactors = FALSE
  ))
}
rows <- list()

for (z in seq_len(nrow(spec))) {
  fam_name <- spec$family[z]
  many <- identical(spec$case[z], "many_visits")
  custom <- identical(spec$case[z], "custom")
  N <- if (custom) as.integer(Sys.getenv("FGEE_BENCH_N")) else 25L
  n_long <- if (custom) {
    as.integer(Sys.getenv("FGEE_BENCH_NI", "100"))
  } else if (many) 100L else 5L
  n_fun <- if (custom) as.integer(Sys.getenv("FGEE_BENCH_L", "50")) else 50L
  p <- if (custom) as.integer(Sys.getenv("FGEE_BENCH_P", "16")) else 16L
  q <- as.integer(Sys.getenv("FGEE_BENCH_Q", "3"))
  fx <- validation_fixture(
    N = N, n_long = n_long, n_fun = n_fun, p = p,
    family = fam_name, seed = 2000L + z
  )
  fit <- validation_fake_fit(p = p, q = q, family = family_object(fam_name))
  exact <- identical(fam_name, "gaussian")
  st <- fastFGEE:::fgee_build_working_stats(
    fx$data, fx$namesd, retain = "scores", exact_gaussian = exact,
    beta0 = fit$coefficients
  )
  prep <- fastFGEE:::fgee_fastk_prepare(
    st, fx$data, fx$namesd, fit, K = 10L, seed = 1L,
    exact = exact, memory = "balanced"
  )
  cv_grid <- list(
    c(1e-2, 1e-1, 1, 10, 100),
    c(0.1, 1, 10),
    c(0.5, 1, 2)
  )

  staged <- validation_time(function() {
    fastFGEE:::fgee_tune_fastk_staged(prep, cv_grid, verbose = FALSE)
  }, repeats)
  qprep <- fastFGEE:::fgee_qreml_prepare(st, fit, exact = exact)
  qreml <- validation_time(function() {
    fastFGEE:::fgee_tune_qreml(qprep, verbose = FALSE)
  }, repeats)
  grad <- validation_time(function() {
    fastFGEE:::fgee_tune_fastk_grad(
      prep, qreml_lambda = qreml$value$lambda,
      start_strategy = "qreml", verbose = FALSE
    )
  }, repeats)

  score_staged <- fastFGEE:::fgee_fastk_score_grad(
    prep, staged$value$lambda, need_gradient = FALSE
  )$score
  score_grad <- fastFGEE:::fgee_fastk_score_grad(
    prep, grad$value$lambda, need_gradient = FALSE
  )$score
  score_qreml <- fastFGEE:::fgee_fastk_score_grad(
    prep, qreml$value$lambda, need_gradient = FALSE
  )$score

  rows[[length(rows) + 1L]] <- data.frame(
    family = fam_name, case = spec$case[z], N = N, n_long = n_long,
    n_fun = n_fun, p = p, M = N * n_long * n_fun,
    staged_sec = staged$median_elapsed,
    gradient_sec = grad$median_elapsed,
    qreml_sec = qreml$median_elapsed,
    staged_over_gradient = staged$median_elapsed / grad$median_elapsed,
    staged_over_qreml = staged$median_elapsed / qreml$median_elapsed,
    staged_score = score_staged,
    gradient_score = score_grad,
    qreml_fastk_score = score_qreml,
    gradient_minus_staged = score_grad - score_staged,
    gradient_evaluations = grad$value$evaluations,
    qreml_evaluations = qreml$value$evaluations,
    phi_pen = qreml$value$information$phi_pen,
    phi_used = qreml$value$phi_used
  )
}

ans <- do.call(rbind, rows)
utils::write.csv(ans, file.path(out_dir, "tuning_benchmark.csv"), row.names = FALSE)
print(ans, row.names = FALSE)
cat("Wrote", file.path(out_dir, "tuning_benchmark.csv"), "\n")
