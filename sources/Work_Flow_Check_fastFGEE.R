library(devtools)
path <- "/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE"
# remove old study strap package
detach("package:fastFGEE", unload=TRUE)
remove.packages("fastFGEE")

setwd(path)

# document and build package
devtools::document()
# devtools::build()

devtools::install(build_vignettes = TRUE)

# Now it will be installed correctly
library(fastFGEE)
vignette("fastFGEE")

#######################
# cleanup sequence after editing
setwd("/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE")
unlink(c("man/fgee.Rd", "NAMESPACE"), force = TRUE)

cleanup_fastFGEE_tree <- function(path = ".") {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)

  file_targets <- c(
    ".DS_Store",
    ".Rhistory",
    ".Rapp.history",
    ".Rbuildignore.txt"
  )

  dir_targets <- c(
    ".git",
    ".Rproj.user",
    "__MACOSX"
  )

  files <- list.files(
    path,
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  dirs <- c(
    file.path(path, dir_targets),
    list.dirs(path, full.names = TRUE, recursive = TRUE)
  )
  dirs <- unique(dirs)

  files_to_delete <- files[basename(files) %in% file_targets]
  dirs_to_delete  <- dirs[basename(dirs) %in% dir_targets]

  if (length(files_to_delete)) {
    unlink(files_to_delete, force = TRUE)
  }
  if (length(dirs_to_delete)) {
    unlink(dirs_to_delete, recursive = TRUE, force = TRUE)
  }

  invisible(list(
    files_deleted = files_to_delete,
    dirs_deleted = dirs_to_delete
  ))
}

# source("cleanup_fastFGEE_tree.R")
cleanup_fastFGEE_tree("/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE")


devtools::document()
devtools::build()
devtools::check()


# install new one
# install.packages("/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE_0.1.0.tar.gz",
# repos = NULL, type = "source")
# library(fastFGEE)

# document and build package -- Need to repeat this again because need smtl_setup run before vignette can be
# properly built
setwd(path)

devtools::document()
devtools::build()

devtools::check()

# STOP HERE !!!!






# clean up suggested by Chat GPT

#============================================
library(devtools)

pkg <- "/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE"

# unload package from current session
if ("package:fastFGEE" %in% search()) {
  detach("package:fastFGEE", unload = TRUE, character.only = TRUE)
}

setwd(pkg)

# At this point, the stale exported functions should already be removed
# or no longer tagged with #' @export:
#   fun.sandwich.ind
#   fun.gee1step.dist
#   fun.gee1step.dist_fullyItr
#   fgee.plot.old

# remove stale generated docs so roxygen rebuilds them cleanly
unlink(file.path(pkg, "NAMESPACE"), force = TRUE)
unlink(file.path(pkg, "man", c(
  "fgee.Rd",
  "fgee.plot.Rd",
  "fun.sandwich.ind.Rd",
  "fun.gee1step.dist.Rd",
  "fun.gee1step.dist_fullyItr.Rd",
  "fgee.plot.old.Rd"
)), force = TRUE)

# regenerate NAMESPACE + man/*.Rd from roxygen comments
devtools::document(pkg = pkg, roclets = c("rd", "namespace"))

# inspect exports
ns_lines <- readLines(file.path(pkg, "NAMESPACE"))
grep("^export\\(", ns_lines, value = TRUE)

#============================================
