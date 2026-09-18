#!/bin/sh
set -eu

SRC_INPUT=${1:-$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)}
SRC=$(CDPATH= cd -- "$SRC_INPUT" && pwd)
OUT_INPUT=${2:-"$(dirname -- "$SRC")/fastFGEE_0.3.0.9006_R_validation"}
OUT_PARENT=$(dirname -- "$OUT_INPUT")
OUT_BASE=$(basename -- "$OUT_INPUT")
mkdir -p "$OUT_PARENT"
OUT_PARENT=$(CDPATH= cd -- "$OUT_PARENT" && pwd)
OUT="$OUT_PARENT/$OUT_BASE"
R_BIN=${R_BIN:-R}
RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}

case "$OUT/" in
  "$SRC/"*)
    echo "ERROR: validation output must be outside the package source tree" >&2
    exit 2
    ;;
esac

if ! command -v "$R_BIN" >/dev/null 2>&1; then
  echo "ERROR: R executable not found: $R_BIN" >&2
  exit 127
fi
if ! command -v "$RSCRIPT_BIN" >/dev/null 2>&1; then
  echo "ERROR: Rscript executable not found: $RSCRIPT_BIN" >&2
  exit 127
fi

rm -rf "$OUT"
mkdir -p "$OUT" "$OUT/library" "$OUT/check"

# Verify that the committed registration files are exactly what the installed
# Rcpp generator produces. Work in a copy so validation never edits the source.
ATTR_COPY="$OUT/attributes_source"
cp -R "$SRC" "$ATTR_COPY"
rm -rf "$ATTR_COPY/.git" "$ATTR_COPY/validation"
"$RSCRIPT_BIN" -e 'args <- commandArgs(TRUE); if (!requireNamespace("Rcpp", quietly=TRUE)) stop("Rcpp required"); Rcpp::compileAttributes(args[[1]])' "$ATTR_COPY"
diff -u "$SRC/R/RcppExports.R" "$ATTR_COPY/R/RcppExports.R" > "$OUT/RcppExports.R.diff" || {
  cat "$OUT/RcppExports.R.diff" >&2
  exit 1
}
diff -u "$SRC/src/RcppExports.cpp" "$ATTR_COPY/src/RcppExports.cpp" > "$OUT/RcppExports.cpp.diff" || {
  cat "$OUT/RcppExports.cpp.diff" >&2
  exit 1
}

# Build from the reconstructed source tree.
(
  cd "$OUT"
  "$R_BIN" CMD build "$SRC"
) > "$OUT/build.log" 2>&1
TARBALL=$(find "$OUT" -maxdepth 1 -type f -name 'fastFGEE_0.3.0.9006.tar.gz' -print | head -n 1)
if [ -z "$TARBALL" ]; then
  echo "ERROR: R CMD build did not create fastFGEE_0.3.0.9006.tar.gz" >&2
  cat "$OUT/build.log" >&2
  exit 1
fi

"$R_BIN" CMD INSTALL --library="$OUT/library" "$TARBALL" > "$OUT/install.log" 2>&1

"$RSCRIPT_BIN" -e 'args <- commandArgs(TRUE); .libPaths(c(args[[1]], .libPaths())); library(fastFGEE); source(args[[2]], echo=TRUE)' \
  "$OUT/library" "$SRC/inst/validation/validate_9006_numerics.R" \
  > "$OUT/focused_validation.log" 2>&1

"$RSCRIPT_BIN" -e 'args <- commandArgs(TRUE); .libPaths(c(args[[1]], .libPaths())); if (!requireNamespace("testthat", quietly=TRUE)) stop("testthat required"); library(fastFGEE); testthat::test_dir(args[[2]], reporter="summary", stop_on_failure=TRUE, stop_on_warning=FALSE)' \
  "$OUT/library" "$SRC/tests/testthat" > "$OUT/testthat.log" 2>&1

(
  cd "$OUT/check"
  _R_CHECK_FORCE_SUGGESTS_=false "$R_BIN" CMD check --as-cran --no-manual "$TARBALL"
) > "$OUT/check.log" 2>&1

CHECK_00=$(find "$OUT/check" -type f -name 00check.log -print | head -n 1)
if [ -n "$CHECK_00" ]; then
  cp "$CHECK_00" "$OUT/00check.log"
  CHECK_STATUS=$(grep '^Status:' "$CHECK_00" | tail -n 1 || true)
else
  CHECK_STATUS="Status line unavailable; inspect check.log"
fi
if printf '%s\n' "$CHECK_STATUS" | grep -q 'ERROR'; then
  echo "ERROR: R CMD check reported an ERROR: $CHECK_STATUS" >&2
  exit 1
fi

{
  printf '%s\n' \
    "Rcpp registration regeneration: PASS" \
    "R CMD build: PASS" \
    "R CMD INSTALL: PASS" \
    "focused numerical validation: PASS" \
    "testthat: PASS" \
    "R CMD check command: COMPLETED"
  printf 'R CMD check summary: %s\n' "${CHECK_STATUS:-OK}"
} > "$OUT/FINAL_STATUS.txt"
cat "$OUT/FINAL_STATUS.txt"
