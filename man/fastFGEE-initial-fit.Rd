% Hand-maintained internal documentation for fastFGEE.
\name{fastFGEE-initial-fit}
\alias{fastFGEE-initial-fit}
\title{Internal baseline inference around a supplied \code{pffr()} fit}
\description{
Notes on the internal initial-fit path used by fastFGEE's development and
validation code. This page documents internal behaviour only. The control that
selects this path is internal and historical; it is not an argument of
\code{\link{fgee}}, and public one-step estimation is unaffected by anything
described here.
}
\section{What the path computes}{
The internal initial-fit path keeps the supplied \code{pffr()} coefficient
estimates fixed. It does not re-estimate the mean. It computes plug-in sandwich
and wild-bootstrap quantities from fastFGEE's own working scores and working
covariance, with the initial fit's smoothing penalty converted to the
cluster-average scale used by the one-step bread \code{Wbar + P}.

mgcv's smoothing parameters penalise the total (unscaled) deviance or residual
sum of squares, so the penalised Hessian is \code{(X'WX + S_lambda) / phi}. The
one-step equations here are on the averaged cluster scale
\code{Wbar = N^-1 sum_i D_i' V_i^-1 D_i}, so the matching penalty is
\code{S_lambda / (N * c_work)}, where \code{N} is the number of independent
clusters and \code{c_work} is the scalar dispersion actually multiplied into the
working variance.
}
\section{Family scale factors}{
\code{c_work} is the dispersion present in the working variance, which is not
always \code{fit$sig2}:
\itemize{
  \item Gaussian: \code{c_work = fit$sig2}. The working variance is the
    empirical residual variance at each functional grid point, whose average is
    close to \code{sig2}; this grid-specific convention is retained.
  \item Gamma: \code{c_work} is the dispersion multiplied into the working
    variance, so the divisor is \code{N} times that dispersion, not \code{N}.
  \item Quasi families on the current fixed-nuisance path: the working variance
    carries no extra scalar dispersion, so the divisor is \code{N} and
    \emph{not} \code{N * fit$sig2}, even when the fitted \code{sig2} is far
    from one.
  \item Binomial and Poisson: \code{c_work = 1}.
  \item Negative binomial and beta: \code{c_work = 1}. The size and precision
    parameters already enter \code{V(mu)} as they do in mgcv's extended
    families, so no further global factor is applied.
}
}
\section{What is and is not rescaled}{
The original \code{pffr()} coefficient estimates are unchanged, and so are the
smoothing parameters stored on the original fitted model object. What this path
\emph{returns} is converted: its \code{lambda} and \code{pen.mat} are on the
cluster-average scale, that is \code{sp / (N * c_work)}. Code comparing those
returned values against \code{fit$sp} must apply the same factor.
}
\section{Scope and limitations}{
The result is not universally identical to the sandwich for the original
\code{pffr()} likelihood estimator, and correct penalty normalisation alone does
not establish estimator-specific covariance validity when the estimating
equations differ:
\itemize{
  \item The Gaussian path uses grid-specific working variances rather than the
    single scalar variance of the original Gaussian \code{pffr()} fit.
  \item For beta regression the penalty scale factor is the number of
    independent clusters, which is correct, but fastFGEE's GEE mean score is not
    the \code{mgcv::betar} likelihood score.
  \item For Gamma the two internal engines may use different working
    dispersions: one takes the dispersion from the nuisance record, while the
    legacy path estimates one from the current fitted mean and warns that it has
    done so. Each engine's returned penalty is consistent with the dispersion it
    used, so comparisons must use that engine's own value.
}
These are properties of this internal baseline path. They are not claims about
the public one-step estimator, whose estimation and smoothing-parameter
selection are unchanged.
}
\keyword{internal}
