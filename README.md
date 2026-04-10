# fastFGEE

`fastFGEE` fits fast one-step functional generalized estimating equations (fGEE)
for longitudinal functional outcomes. The package uses a `refund::pffr()` initial
fit and then updates the coefficient estimate with a working correlation structure
in the longitudinal and/or functional direction.

## Main features

- one-step penalized fGEE estimation
- supports quasi-likelihoods derived from many families (e.g., Gaussian, binomial, Poisson, Gamma, negative-binomial, beta) through the package's family-handling utilities
- supports standard link functions such as identity, log, logit, probit, inverse
- supports working covariances that can model working correlations in longitudinal and/or functional directions as independent, exchangeable, AR1, and FPCA-based 
- fast cluster cross-validation for smoothing-parameter selection
- sandwich and bootstrap-based uncertainty quantification 
- pointwise and joint confidence intervals that yield valid inference even when the working correlation is misspecified


## Installation

### Github version

```r
install.packages("fastFGEE")
```
### Development version

```r
remotes::install_github("gloewing/fastFGEE")
```


### Optional accelerator

The package works without `sanic`, but matrix solves can be faster when it is
installed:

```r
install.packages("sanic")
```

### Optional archived package for irregular AR(1) work

If you choose to support archived `irregulAR1`-based workflows outside CRAN,
install it manually from the CRAN archive:

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/irregulAR1/irregulAR1_1.0.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

## Example

```r
library(fastFGEE)
data("DTI", package = "refund")

set.seed(1)
ids <- sample(unique(DTI$ID), 10)
DTI_use <- subset(DTI, ID %in% ids)
DTI_use <- data.frame(
  cca = I(DTI_use$cca),
  case = DTI_use$case,
  visit = as.numeric(DTI_use$visit),
  sex = DTI_use$sex,
  ID = DTI_use$ID
)

fit <- fgee(
  formula = cca ~ case + visit,
  data = DTI_use,
  cluster = "ID",
  family = gaussian(),
  corr_fn = "independent",
  corr_long = "exchangeable"
)

fgee.plot(fit)
```

## Development note

Portions of the package were developed with assistance from large language
models. All code was reviewed and validated by the package author.
