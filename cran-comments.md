## Test environments

* local macOS: Apple M1 Max, macOS 15.7.2, R 4.5.1
* win-builder (release): R 4.5.3
* win-builder (devel): R-devel (2026-03-25 r89703)

## R CMD check results

0 errors | 0 warnings | 0 notes on local `devtools::check()` / `R CMD check --as-cran`

win-builder checks on release and devel completed with no package warnings or errors.
The only remaining message is the CRAN incoming NOTE concerning the suggested package `irregulAR1`.

## Notes

This is a new submission.

`irregulAR1` is listed only in `Suggests` and is used conditionally for
irregularly spaced AR(1) precision matrices. The package does not require
`irregulAR1` for installation, loading, examples, tests, or vignettes.
Standard functionality works without it.

As recommended by CRAN policy, the `Description` field explains that this is
an optional archived package obtained from the CRAN Archive.

Package URLs were rechecked before submission.