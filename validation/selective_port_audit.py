#!/usr/bin/env python3
"""Verify that 0.3.0.9006 is a bounded port from the validated 0.3.0.9004 tree."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "df79730"

ALLOWED_EXACT = {
    ".Rbuildignore",
    "DESCRIPTION",
    "DEVELOPMENT.md",
    "NAMESPACE",
    "NEWS.md",
    "R/RcppExports.R",
    "R/WD_estimate.R",
    "R/cv.R",
    "R/family_fns.R",
    "R/fastFGEE-package.R",
    "R/fgee.R",
    "R/fgee.plot.R",
    "R/fgee_public.R",
    "R/irregular_ar1_internal.R",
    "R/linear_algebra_internal.R",
    "R/nuisance_integration_core.R",
    "R/solve_pd.R",
    "R/working_stats.R",
    "R/zzz_nuisance_integration.R",
    "cran-comments.md",
    "inst/validation/README.md",
    "inst/validation/benchmark_fgee_engines.R",
    "inst/validation/validate_9006_numerics.R",
    "man/fastFGEE-package.Rd",
    "man/fgee.Rd",
    "man/fgee.plot.Rd",
    "src/Makevars",
    "src/Makevars.win",
    "src/RcppExports.cpp",
    "src/fgee_linear_algebra.cpp",
    "tests/testthat/test-dependency-and-registration.R",
    "tests/testthat/test-fgee-integration.R",
    "tests/testthat/test-fgee-newfamilies.R",
    "tests/testthat/test-internal-linear-algebra.R",
    "tests/testthat/test-public-estimator-surface.R",
    "vignettes/fastFGEE.Rmd",
}
ALLOWED_PREFIXES = ("validation/",)
PRESERVED = (
    "R/tuning_fastk.R",
    "R/tuning_dispatch.R",
    "R/nuisance_parameters.R",
)


def run(*args: str) -> str:
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    failures: list[str] = []
    # Include both tracked differences and untracked development files.
    changed = set(filter(None, run("git", "diff", "--name-only", BASE, "--").splitlines()))
    # Do not strip porcelain output: the leading status-space on a worktree-only
    # modification is significant (for example, " M .Rbuildignore").
    status = subprocess.check_output(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=ROOT,
        text=True,
    )
    for line in status.splitlines():
        path = line[3:]
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        changed.add(path)

    unexpected = sorted(
        path for path in changed
        if path not in ALLOWED_EXACT and not path.startswith(ALLOWED_PREFIXES)
    )
    if unexpected:
        failures.append(f"unexpected changed paths: {unexpected}")

    preserved: dict[str, dict[str, object]] = {}
    for path in PRESERVED:
        base_bytes = subprocess.check_output(
            ["git", "show", f"{BASE}:{path}"], cwd=ROOT
        )
        current_bytes = (ROOT / path).read_bytes()
        same = base_bytes == current_bytes
        preserved[path] = {
            "baseline_sha256": sha256(base_bytes),
            "current_sha256": sha256(current_bytes),
            "byte_identical": same,
        }
        if not same:
            failures.append(f"validated 0.3.0.9004 file changed unexpectedly: {path}")

    required_new = {
        "R/fgee_public.R",
        "R/irregular_ar1_internal.R",
        "R/linear_algebra_internal.R",
        "R/nuisance_integration_core.R",
        "src/fgee_linear_algebra.cpp",
        "inst/validation/validate_9006_numerics.R",
    }
    missing_new = sorted(path for path in required_new if not (ROOT / path).is_file())
    if missing_new:
        failures.append(f"required selective-port files are missing: {missing_new}")
    if (ROOT / "R/zzz_nuisance_integration.R").exists():
        failures.append("load-order nuisance wrapper was not removed")

    result = {
        "baseline_commit": BASE,
        "changed_paths": sorted(changed),
        "unexpected_changed_paths": unexpected,
        "preserved_9004_files": preserved,
        "required_new_files_missing": missing_new,
        "failures": failures,
        "passed": not failures,
    }
    out = ROOT / "validation/selective_port_audit.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if not failures else 1)


if __name__ == "__main__":
    main()
