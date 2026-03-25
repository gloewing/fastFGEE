#' fastFGEE: Fast Functional Generalized Estimating Equations
#'
#' Fits functional generalized estimating equations (fGEE) for
#' longitudinal functional outcomes using a one-step estimator.
#'
#' The main user-facing functions are \code{\link{fgee}} and \code{\link{fgee.plot}}.
#'
#' @section Optional packages:
#' \itemize{
#'   \item \code{sanic} is used when available for faster positive-definite solves.
#'   \item \code{irregulAR1} is used only for irregularly spaced AR1 precision matrices.
#'   \item \code{Rcpp} and \code{RcppArmadillo} are used only for developer-side experimental
#'   C++ loss code and are not required for the minimal CRAN build.
#' }
#'
#' @section AI-assisted development:
#' Portions of the package were developed with assistance from large language
#' models. All code was reviewed, tested, and edited by the human author(s),
#' who take responsibility for the package.
#'
#' @keywords internal
#' @import data.table
#' @importFrom ggplot2 aes coord_cartesian element_text geom_hline
#'   geom_line geom_ribbon geom_segment ggplot labs
#'   scale_colour_manual theme theme_classic
#' @importFrom mgcv s
#' @importFrom stats na.omit
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
