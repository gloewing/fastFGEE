#' fastFGEE: Fast Functional Generalized Estimating Equations
#'
#' Fits functional generalized estimating equations (fGEE) for longitudinal
#' functional outcomes using the validated one-step estimator.
#'
#' The main user-facing functions are [fgee()] and [fgee.plot()].
#'
#' @section Supported families:
#' The package supports quasi-likelihoods derived from a range of
#' distributions. Substantial simulation evidence exists for quasi-likelihoods
#' derived from Gaussian, binomial, Poisson, negative binomial, Gamma and beta
#' families, across a range of sample sizes, cluster sizes and working
#' correlation specifications. Other families understood by `refund::pffr()`
#' are likely to be supported but have not been examined as thoroughly in
#' simulation, so results for them should be interpreted with corresponding
#' care. For negative binomial outcomes use `mgcv::nb()`; for proportion
#' outcomes use `mgcv::betar()`.
#'
#' @section Numerical backends:
#' The package directly imports Rcpp and links its registered compiled routines.
#' Internal LAPACK Cholesky routines are used for symmetric positive-definite
#' solves. Exact tridiagonal Markov precision operations for irregularly sampled
#' continuous-time AR(1) working correlations are implemented within fastFGEE
#' from the formulas of Allevius (2018); no code from the archived
#' `irregulAR1` package is used. `SuperGauss` remains an optional backend for
#' explicitly requested regular-grid Toeplitz calculations. `RcppArmadillo` is
#' retained only for historical developer-side experimental CV helpers.
#'
#' @section Inference scope:
#' Pointwise intervals and optional simultaneous bands answer different
#' inferential questions. A requested simultaneous band uses a whole-curve
#' maximum statistic, but nominal finite-sample simultaneous coverage is not
#' guaranteed, particularly with few or high-leverage independent clusters.
#'
#' @section AI-assisted development:
#' Portions of the package were developed with assistance from large language
#' models. All code remains the responsibility of the human author(s).
#'
#' @references
#' Allevius, B. (2018). On the precision matrix of an irregularly sampled AR(1)
#' process. *Stockholm University Research Report*.
#'
#' Loewinger, G., Levis, A. W., Cui, E., and Pereira, F. (2025). Fast Penalized
#' Generalized Estimating Equations for Large Longitudinal Functional Datasets.
#'
#' @keywords internal
#' @import data.table
#' @importFrom ggplot2 aes coord_cartesian element_text geom_hline
#'   geom_line geom_ribbon geom_segment ggplot labs
#'   scale_colour_manual theme theme_classic
#' @importFrom mgcv s
#' @importFrom Rcpp evalCpp
#' @importFrom stats na.omit
#' @importFrom stats setNames
#' @importFrom utils head tail
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
