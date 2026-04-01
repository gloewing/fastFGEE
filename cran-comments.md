## Resubmission

This is a resubmission. I addressed the issues from the previous CRAN review as follows:

* added method references directly to the `Description` field in the requested
  format using DOI / URL links;
* made the optional archived package `irregulAR1` access path explicit in the
  `Description` field via the CRAN Archive URL;
* replaced `\dontrun{}` with `\donttest{}` in `fgee` examples;
* removed unsuppressible console output from `R/WD_estimate.R` by replacing
  `print()` calls with informative `stop()` messages.

## Test environments

* local macOS: Apple M1 Max, macOS 15.7.2, R 4.5.1
* win-builder (release): R 4.5.3
* win-builder (devel): R-devel (2026-03-31 r89747)

## R CMD check results

0 errors | 0 warnings | 0 notes on local `devtools::check()` / `R CMD check --as-cran`

win-builder checks on release and devel completed with no package warnings or errors.
The only remaining messages are the CRAN incoming NOTE concerning the suggested
package `irregulAR1` and the automatic DESCRIPTION spell-check note for the
bibliographic reference `Loewinger et al. (2025)`.

## Notes

`irregulAR1` is listed only in `Suggests` and is used conditionally for
optional irregularly spaced AR(1) precision-matrix code paths. The package does
not require `irregulAR1` for installation, loading, examples, tests, or
vignettes. Core functionality works without it.

The DESCRIPTION spell-check NOTE flags `Loewinger`, `et`, and `al` from the
author-year citation `Loewinger et al. (2025)` included in the Description
field to document the method reference requested in the previous review.
This is a bibliographic citation, not a misspelling.