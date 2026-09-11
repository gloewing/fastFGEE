.source_root <- function() {
  test_dir <- normalizePath(testthat::test_path(), mustWork = TRUE)
  candidates <- c(
    file.path(test_dir, "..", ".."),
    file.path(test_dir, "..", "..", "00_pkg_src", "fastFGEE"),
    file.path(getwd(), "00_pkg_src", "fastFGEE"),
    getwd()
  )
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  for (candidate in candidates) {
    desc <- file.path(candidate, "DESCRIPTION")
    if (file.exists(desc) && dir.exists(file.path(candidate, "R")) &&
        dir.exists(file.path(candidate, "src"))) {
      first <- read.dcf(desc)[1, "Package"]
      if (identical(unname(first), "fastFGEE")) return(candidate)
    }
  }
  stop("Could not locate the fastFGEE source tree for source-level tests.")
}

test_that("archived numerical dependencies are absent from executable source", {
  root <- .source_root()
  desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
  dependency_text <- paste(
    desc[intersect(c("Depends", "Imports", "Suggests", "LinkingTo"), names(desc))],
    collapse = " "
  )
  expect_false(grepl("irregulAR1|sanic", dependency_text))
  expect_match(desc[["Imports"]], "Rcpp")
  expect_match(desc[["LinkingTo"]], "Rcpp")
  expect_match(desc[["Imports"]], "refund \\(>= 0\\.1-40\\)")

  source_files <- c(
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
    list.files(file.path(root, "src"), pattern = "[.](cpp|h)$", full.names = TRUE)
  )
  source_text <- paste(unlist(lapply(source_files, readLines, warn = FALSE),
                              use.names = FALSE), collapse = "\n")
  expect_false(grepl("irregulAR1::|sanic::", source_text, fixed = FALSE))
  expect_false(grepl("requireNamespace\\(\"(irregulAR1|sanic)\"", source_text))

  r_source <- paste(unlist(lapply(
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
    readLines, warn = FALSE
  ), use.names = FALSE), collapse = "\n")
  expect_false(grepl(
    "if[[:space:]]*\\([[:space:]]*!?TRUE[[:space:]]*\\)",
    r_source
  ))
})

test_that("every Rcpp attribute export is registered", {
  root <- .source_root()
  cpp <- list.files(file.path(root, "src"), pattern = "[.]cpp$", full.names = TRUE)
  source_text <- paste(unlist(lapply(cpp, readLines, warn = FALSE),
                              use.names = FALSE), collapse = "\n")
  # Source-level expected set is deliberately explicit so this test cannot pass
  # merely because both sides of an automated count omitted the same routine.
  expected <- c(
    "_fastFGEE_fgee_kron_inverse_kernel" = 7L,
    "_fastFGEE_fastk_fold_kernel" = 9L,
    "_fastFGEE_fgee_sympd_inverse_cpp" = 1L,
    "_fastFGEE_fgee_sympd_solve_cpp" = 2L,
    "_fastFGEE_fgee_iar1_precision_bands_cpp" = 2L,
    "_fastFGEE_fgee_iar1_precision_cpp" = 2L,
    "_fastFGEE_fgee_iar1_apply_precision_cpp" = 3L,
    "_fastFGEE_fgee_iar1_profile_nll_cpp" = 3L
  )
  exports <- readLines(file.path(root, "R", "RcppExports.R"), warn = FALSE)
  registrations <- readLines(file.path(root, "src", "RcppExports.cpp"),
                             warn = FALSE)
  for (symbol in names(expected)) {
    expect_true(any(grepl(symbol, exports, fixed = TRUE)), info = symbol)
    expect_true(any(grepl(paste0('"', symbol, '"'), registrations,
                          fixed = TRUE)), info = symbol)
  }
  registration_lines <- grep(
    '^    \\{"_fastFGEE_', registrations, value = TRUE
  )
  expect_length(registration_lines, length(expected))
  for (symbol in names(expected)) {
    line <- registration_lines[grepl(paste0('"', symbol, '"'),
                                     registration_lines, fixed = TRUE)]
    # expect_length() takes no `info`; expect_equal() does, and keeping it
    # is what names the offending symbol when this fails.
    expect_equal(length(line), 1L, info = symbol)
    observed <- as.integer(sub('.*,[[:space:]]*([0-9]+)[[:space:]]*\\},.*',
                               '\\1', line))
    expect_identical(observed, unname(expected[[symbol]]), info = symbol)
  }
  expect_true(grepl("fgee_sympd_solve_cpp", source_text, fixed = TRUE))
})

test_that("nuisance integration has one active definition per core wrapper", {
  root <- .source_root()
  r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
  lines <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
  for (name in c("fgee_update_working_cols_dt", "fgee_build_working_stats",
                 "get_family_info", "fgee")) {
    pattern <- paste0("^", name, "[[:space:]]*<-[[:space:]]*function")
    expect_equal(sum(grepl(pattern, lines)), 1L, info = name)
  }
  expect_false(file.exists(file.path(root, "R", "zzz_nuisance_integration.R")))
  expect_true(is.function(fastFGEE:::.fgee_update_working_cols_dt_core))
  expect_true(is.function(fastFGEE:::.fgee_build_working_stats_core))
  expect_true(is.function(fastFGEE:::.fgee_get_family_info_core))
})
