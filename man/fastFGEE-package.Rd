% Generated from R/fastFGEE-package.R for fastFGEE 0.3.0.9006
\docType{package}
\name{fastFGEE-package}
\alias{fastFGEE}
\alias{fastFGEE-package}
\title{fastFGEE: Fast Functional Generalized Estimating Equations}
\description{
Fits functional generalized estimating equations for longitudinal functional outcomes using the validated one-step estimator.
}
\details{
The main user-facing functions are \code{\link{fgee}} and \code{\link{fgee.plot}}.
}
\section{Supported families}{
The package supports quasi-likelihoods derived from a range of distributions. Substantial simulation evidence exists for quasi-likelihoods derived from Gaussian, binomial, Poisson, negative binomial, Gamma and beta families, across a range of sample sizes, cluster sizes and working correlation specifications. Other families understood by \code{refund::pffr()} are likely to be supported but have not been examined as thoroughly in simulation, so results for them should be interpreted with corresponding care. For negative binomial outcomes use \code{mgcv::nb()}; for proportion outcomes use \code{mgcv::betar()}.
}

\section{Numerical backends}{
The package directly imports Rcpp and links registered compiled routines. Internal LAPACK Cholesky routines are used for symmetric positive-definite solves. Exact tridiagonal Markov precision operations for irregularly sampled continuous-time AR(1) working correlations are implemented within fastFGEE from the formulas of Allevius (2018); no code from the archived \pkg{irregulAR1} package is used. \pkg{SuperGauss} remains an optional backend for explicitly requested regular-grid Toeplitz calculations. \pkg{RcppArmadillo} is retained only for historical developer-side experimental CV helpers.
}
\section{Inference scope}{
Pointwise intervals and optional simultaneous bands answer different inferential questions. A requested simultaneous band uses a whole-curve maximum statistic, but nominal finite-sample simultaneous coverage is not guaranteed, particularly with few or high-leverage independent clusters.
}
\section{AI-assisted development}{
Portions of the package were developed with assistance from large language models. All code remains the responsibility of the human author(s).
}
\references{
Allevius, B. (2018). On the precision matrix of an irregularly sampled AR(1) process.

Loewinger, G., Levis, A. W., Cui, E., and Pereira, F. (2025). Fast Penalized Generalized Estimating Equations for Large Longitudinal Functional Datasets.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/gloewing/fastFGEE}
  \item Report bugs at \url{https://github.com/gloewing/fastFGEE/issues}
}
}
\author{
\strong{Maintainer}: Gabriel Loewinger \email{gloewinger@gmail.com} (\href{https://orcid.org/0000-0002-0755-8520}{ORCID})
}
\keyword{internal}
