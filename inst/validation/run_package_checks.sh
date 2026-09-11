#!/usr/bin/env bash
set -euo pipefail

PKG_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
OUT=${FGEE_CHECK_OUT:-"$PKG_ROOT/check-output"}
LIB=${R_LIBS_USER:-"$OUT/Rlib"}
mkdir -p "$OUT" "$LIB"
export R_LIBS_USER="$LIB"

cd "$(dirname "$PKG_ROOT")"
R CMD build --no-build-vignettes "$PKG_ROOT"
TARBALL=$(ls -t fastFGEE_*.tar.gz | head -1)
R CMD check --no-manual --no-build-vignettes \
  --library="$LIB" --output="$OUT" "$TARBALL"

echo "Built: $(pwd)/$TARBALL"
echo "Check output: $OUT"
