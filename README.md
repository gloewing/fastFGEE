# fastFGEE

`fastFGEE` fits fast one-step functional generalized estimating equations (fGEE)
for longitudinal functional outcomes. The package uses a `refund::pffr()` initial
fit and then updates the coefficient estimate with a working correlation structure
in the longitudinal and/or functional direction. See the [vignette](https://cran.r-project.org/web/packages/fastFGEE/vignettes/fastFGEE.html) for examples.

## Main features

- one-step penalized fGEE estimation
- supports quasi-likelihoods derived from many families (e.g., Gaussian, binomial, Poisson, Gamma, negative-binomial, beta) through the package's family-handling utilities
- supports standard link functions such as identity, log, logit, probit, inverse
- supports working covariances that can model working correlations in longitudinal and/or functional directions as independent, exchangeable, AR1, and FPCA-based 
- fast cluster cross-validation for smoothing-parameter selection
- sandwich and bootstrap-based uncertainty quantification 
- pointwise and joint confidence intervals that yield valid inference even when the working correlation is misspecified


## Installation

### From r-universe (recommended -- no compiler needed)

```r
install.packages("fastFGEE", repos = c(
  gloewing = "https://gloewing.r-universe.dev",
  CRAN     = "https://cloud.r-project.org"
))
```

This installs a prebuilt binary on Windows and macOS, so no development tools
are required.

### CRAN version

```r
install.packages("fastFGEE")
```

### Development version from source

`fastFGEE` contains compiled C++ code, so installing from source requires a
toolchain: **Rtools** on Windows, **Xcode Command Line Tools** on macOS
(`xcode-select --install`), or `build-essential` on Linux.

```r
remotes::install_github("gloewing/fastFGEE")
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
