% Generated from R/fgee_public.R for fastFGEE 0.3.0.9006
\name{fgee}
\alias{fgee}
\title{Fit a One-Step Functional Generalized Estimating Equation}
\usage{
fgee(
  formula,
  data,
  cluster,
  family,
  corr_fn = "ar1",
  corr_long = "ar1",
  time = NULL,
  long.dir = TRUE,
  var.type = "sandwich",
  pffr.mod = NULL,
  knots = NULL,
  bs = "bs",
  cv = "fastkfold",
  cv.grid = NULL,
  rho.smooth = FALSE,
  rho.pool = c("fn", "none", "long", "both"),
  joint.CI = "wild",
  linpred_method = c("accumulate", "matrix"),
  clip_mu = 0,
  m.pffr = c(2, 1),
  check_alignment = TRUE,
  boot.samps = 3000,
  sp.method = c("auto", "fastk_staged", "fastk_grad", "fastk_grad_fast",
    "sandwich_qreml", "qreml_fastk"),
  working.retain = c("auto", "scores", "aggregate", "full"),
  corr.solver = c("auto", "exact", "supergauss"),
  fastk.K = 10L,
  fastk.seed = 1L,
  fastk.memory = c("balanced", "speed", "lowmem"),
  fastk.start = c("qreml", "robust", "compact"),
  fastk.kernel = fgee_fastk_kernel_ok(),
  qreml.phi.method = c("penalized", "all", "fixed"),
  qreml.phi.fixed = NULL,
  qreml.phi.weight = 1,
  qreml.phi.clip = c(1, 8),
  keep.tuning.workspace = FALSE,
  verbose.tuning = TRUE,
  keep.data = TRUE,
  keep.initial.fit = TRUE,
  keep.working.stats = TRUE,
  ...
)
}
\arguments{
\item{formula}{A model formula whose left-hand side is a functional response.}
\item{data}{A data frame containing the variables in \code{formula}.}
\item{cluster}{Name of the cluster identifier column.}
\item{family}{A family object or family name understood by \code{refund::pffr()}. Substantial simulation evidence exists for quasi-likelihoods derived from Gaussian, binomial, Poisson, negative binomial, Gamma and beta families; other families supported by \code{refund::pffr()} are likely to work but have not been examined as thoroughly. For negative binomial outcomes use \code{mgcv::nb()}, and for proportion outcomes use \code{mgcv::betar()}; \code{MASS::negative.binomial()} is not accepted and raises an informative error.}
\item{corr_fn}{Working correlation in the functional direction: \code{"independent"}, \code{"exchangeable"}, \code{"ar1"}, or \code{"fpca"}.}
\item{corr_long}{Working correlation in the longitudinal direction: \code{"independent"}, \code{"exchangeable"}, or \code{"ar1"}.}
\item{time}{Optional name of the longitudinal ordering variable.}
\item{long.dir}{Logical retained for backwards compatibility.}
\item{var.type}{Variance estimator: \code{"sandwich"}, \code{"fastboot"}, or \code{"boot"}.}
\item{pffr.mod}{Optional fitted \code{refund::pffr()} object used as the initial estimator.}
\item{knots}{Number of spline knots for the initial \code{pffr()} fit.}
\item{bs}{Basis type passed to \code{refund::pffr()}.}
\item{cv}{Cross-validation mode used for smoothing-parameter selection. The public one-step interface currently supports only \code{"fastkfold"}.}
\item{cv.grid}{Optional grid or staged grids of smoothing parameters.}
\item{rho.smooth}{Logical; smooth pointwise working-correlation estimates.}

\item{rho.pool}{Controls pooling of the working-correlation parameter when only one direction is correlated. \code{"fn"} (default) estimates a single pooled functional correlation when \code{corr_long = "independent"}; \code{"long"} pools the longitudinal correlation when \code{corr_fn = "independent"}; \code{"both"} pools either; \code{"none"} lets the surviving parameter vary over the opposite index, reproducing the behaviour of earlier versions. When both directions are correlated the working correlation is separable and both parameters are pooled regardless of this setting.}
\item{joint.CI}{Controls interval construction. \code{"wild"} requests the studentized wild-cluster procedure and its optional simultaneous band. A requested band is a simultaneous-inference procedure, not a guarantee of nominal finite-sample coverage.}
\item{linpred_method}{Method used to form linear predictors.}
\item{clip_mu}{Lower bound used for numerical stabilization of fitted means.}
\item{m.pffr}{Penalty-order specification passed to \code{refund::pffr()}.}
\item{check_alignment}{Logical; check alignment between the wide data and the long representation produced by \code{pffr()}.}
\item{boot.samps}{Number of bootstrap replicates when applicable.}
\item{sp.method}{Smoothing-parameter selector. \code{"auto"} uses staged fastK for Gaussian identity-link models and \code{"fastk_grad_fast"} otherwise. Explicit choices are \code{"fastk_staged"}, \code{"fastk_grad"}, \code{"fastk_grad_fast"}, \code{"sandwich_qreml"}, and \code{"qreml_fastk"}. The last is retained for reproducibility but normally reaches the same exact fastK solution by a slower route.}
\item{working.retain}{Retention profile for compact working statistics: \code{"auto"}, \code{"scores"}, \code{"aggregate"}, or \code{"full"}.}
\item{corr.solver}{Correlation inverse backend. \code{"auto"} uses direct package operators; \code{"supergauss"} forces the optional Toeplitz backend.}
\item{fastk.K}{Number of cluster folds for optimized fastK tuning.}
\item{fastk.seed}{Seed used to construct optimized fastK folds.}
\item{fastk.memory}{Non-Gaussian fastK workspace. \code{"balanced"} is recommended and is the only layout eligible for the compiled loss/gradient kernel. \code{"lowmem"} is deprecated because validation found no peak-memory advantage and it cannot use the compiled kernel.}
\item{fastk.start}{Initialization strategy for analytic-gradient fastK.}
\item{fastk.kernel}{Logical; allow the compiled loss/gradient kernel when supported.}
\item{qreml.phi.method}{Information-scaling diagnostic used by sandwich qREML.}
\item{qreml.phi.fixed}{Fixed information scale when \code{qreml.phi.method = "fixed"}.}
\item{qreml.phi.weight}{Weight applied to the information mismatch before clipping.}
\item{qreml.phi.clip}{Lower and upper safeguards for the qREML information scale.}
\item{keep.tuning.workspace}{Logical; retain large tuning workspaces for debugging.}
\item{verbose.tuning}{Logical; print optimized tuning progress.}
\item{keep.data}{Logical; retain long-format working data in the fitted object.}
\item{keep.initial.fit}{Logical; retain the initial \code{pffr()} fit.}
\item{keep.working.stats}{Logical; retain compact initial and final working statistics.}
\item{...}{Reserved to detect removed estimator controls. Unknown entries are rejected rather than silently ignored.}
}
\value{
An object of class \code{"fgee1step"}. Important components include \code{beta}, \code{vb}, \code{model}, \code{lambda}, and \code{tuning}.
}
\description{
Fits longitudinal functional-response generalized estimating equations using the validated one-step estimator.
}
\details{
Public calls always use one coefficient update and the optimized working-statistics engine. FastK selectors use the one-step fast cluster cross-validation construction; \code{"sandwich_qreml"} remains an explicit opt-in selector rather than a cross-validation method. Historical exact-GLS, legacy-engine, pffr-only, and fully iterated controls remain internal for regression testing and method development and are not accepted by \code{fgee()}.

Pointwise intervals and simultaneous bands answer different inferential questions. A simultaneous band uses a whole-curve maximum statistic in the asymptotic or resampling procedure. Nominal finite-sample simultaneous coverage is not guaranteed, particularly with few or high-leverage independent clusters or unstable nuisance and correlation estimates.
}
\examples{
\dontrun{
data("d", package = "fastFGEE")
fit <- fgee(
  Y ~ X1 + X2,
  data = d,
  cluster = "ID",
  family = binomial(link = "logit"),
  time = "time",
  corr_long = "exchangeable",
  corr_fn = "independent"
)
fgee.plot(fit)
}
}
\references{
Loewinger, G., Levis, A. W., Cui, E., and Pereira, F. (2025). Fast Penalized Generalized Estimating Equations for Large Longitudinal Functional Datasets.
}
\seealso{\code{\link{fgee.plot}}}
