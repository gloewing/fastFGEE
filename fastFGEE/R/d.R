#' Simulated longitudinal functional example data
#'
#' A simulated dataset used in examples and testing for \pkg{fastFGEE}.
#' The object \code{d} is a data frame containing a functional response stored
#' in an \code{AsIs} matrix column together with a cluster identifier, two
#' scalar covariates, and a longitudinal time variable.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'   \item{\code{Y}}{An \code{AsIs} matrix-valued column containing the
#'   functional response. In the included example, the matrix has 100 columns
#'   named \code{Y_1} to \code{Y_100}.}
#'   \item{\code{ID}}{Cluster identifier.}
#'   \item{\code{X1}}{First scalar covariate.}
#'   \item{\code{X2}}{Second scalar covariate.}
#'   \item{\code{time}}{Longitudinal time variable.}
#' }
#'
#' @details
#' This dataset is intended for package examples, vignettes, and quick testing
#' of \code{\link{fgee}}. It represents a simulated binary-response setting on
#' a common functional grid.
#'
#' @usage data(d)
#' @docType data
#' @name d
#'
#' @source Simulated for the package examples.
"d"
