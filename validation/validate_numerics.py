#!/usr/bin/env python3
"""Independent numerical validation for fastFGEE 0.3.0.9006 kernels.

This does not compile or execute Rcpp. It independently checks the formulas
ported into src/fgee_linear_algebra.cpp against NumPy dense linear algebra.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


def bands(time: np.ndarray, rho: float) -> tuple[np.ndarray, np.ndarray]:
    time = np.asarray(time, dtype=float)
    if time.ndim != 1 or len(time) < 1 or np.any(~np.isfinite(time)):
        raise ValueError("bad time")
    if len(time) > 1 and np.any(np.diff(time) <= 0):
        raise ValueError("time not strictly increasing")
    if not (0 <= rho < 1):
        raise ValueError("bad rho")
    n = len(time)
    diag = np.zeros(n)
    off = np.zeros(max(0, n - 1))
    if n == 1:
        diag[0] = 1.0
        return diag, off
    if rho == 0:
        a = np.zeros(n - 1)
        den = np.ones(n - 1)
    else:
        loga = math.log(rho) * np.diff(time)
        a = np.exp(loga)
        den = -np.expm1(2 * loga)
    diag[0] = 1 / den[0]
    if n > 2:
        diag[1:-1] = 1 / den[:-1] + a[1:] ** 2 / den[1:]
    diag[-1] = 1 / den[-1]
    off[:] = -a / den
    return diag, off


def precision(time: np.ndarray, rho: float) -> np.ndarray:
    d, o = bands(time, rho)
    q = np.diag(d)
    if len(o):
        q += np.diag(o, 1) + np.diag(o, -1)
    return q


def apply_precision(rhs: np.ndarray, time: np.ndarray, rho: float) -> np.ndarray:
    x = np.asarray(rhs, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    d, o = bands(time, rho)
    out = d[:, None] * x
    if len(o):
        out[:-1] += o[:, None] * x[1:]
        out[1:] += o[:, None] * x[:-1]
    return out


def profile_nll(resid: np.ndarray, time: np.ndarray, rho: float) -> float:
    resid = np.asarray(resid, dtype=float)
    if rho == 0:
        a = np.zeros(len(time) - 1)
        den = np.ones(len(time) - 1)
    else:
        loga = math.log(rho) * np.diff(time)
        a = np.exp(loga)
        den = -np.expm1(2 * loga)
    q = resid[0] ** 2 + np.sum((resid[1:] - a * resid[:-1]) ** 2 / den)
    return len(resid) * math.log(q / len(resid)) + float(np.sum(np.log(den)))


def main() -> None:
    rng = np.random.default_rng(9006)
    checks: dict[str, float | int | bool] = {}

    # SPD inverse/solve sequence represented by the same mathematical result.
    max_inv = 0.0
    max_solve = 0.0
    for n in (1, 2, 5, 20):
        z = rng.normal(size=(n + 4, n))
        a = z.T @ z + np.eye(n) * 0.25
        b = rng.normal(size=(n, 4))
        inv = np.linalg.inv(a)
        sol = np.linalg.solve(a, b)
        max_inv = max(max_inv, float(np.max(np.abs(a @ inv - np.eye(n)))))
        max_solve = max(max_solve, float(np.max(np.abs(a @ sol - b))))
    checks["spd_inverse_identity_max_abs"] = max_inv
    checks["spd_solve_residual_max_abs"] = max_solve

    max_q = 0.0
    max_apply = 0.0
    max_obj = 0.0
    cases = 0
    for n in (1, 2, 3, 8, 25):
        for rho in (1e-4, 0.12, 0.5, 0.91, 0.999):
            time = np.cumsum(np.r_[0.0, rng.uniform(0.02, 2.5, max(0, n - 1))])
            r = rho ** np.abs(time[:, None] - time[None, :])
            q = precision(time, rho)
            q_ref = np.linalg.inv(r)
            max_q = max(max_q, float(np.max(np.abs(q - q_ref))))
            rhs = rng.normal(size=(n, 5))
            max_apply = max(max_apply, float(np.max(np.abs(apply_precision(rhs, time, rho) - q_ref @ rhs))))
            if n >= 2:
                e = rng.normal(size=n)
                dense_q = float(e @ q_ref @ e)
                sign, logdet = np.linalg.slogdet(r)
                assert sign > 0
                dense_obj = n * math.log(dense_q / n) + logdet
                max_obj = max(max_obj, abs(profile_nll(e, time, rho) - dense_obj))
            cases += 1
    checks["iar1_precision_max_abs"] = max_q
    checks["iar1_apply_max_abs"] = max_apply
    checks["iar1_profile_objective_max_abs"] = max_obj
    checks["iar1_cases"] = cases

    # The zero-correlation boundary is exact working independence. Test it
    # separately so adding the boundary check does not perturb the established
    # random dense-reference cases above.
    zero_time = np.array([0.0, 0.1, 0.8, 2.5])
    zero_rhs = np.arange(8.0).reshape(4, 2)
    checks["iar1_zero_precision_max_abs"] = float(np.max(
        np.abs(precision(zero_time, 0.0) - np.eye(len(zero_time)))
    ))
    checks["iar1_zero_apply_max_abs"] = float(np.max(
        np.abs(apply_precision(zero_rhs, zero_time, 0.0) - zero_rhs)
    ))

    # Exact reduction to regular AR(1) precision on unit gaps.
    max_regular = 0.0
    for n in (2, 3, 20):
        for rho in (0.1, 0.6, 0.95):
            q = precision(np.arange(n, dtype=float), rho)
            ref = np.zeros((n, n))
            den = 1 - rho**2
            ref[0, 0] = ref[-1, -1] = 1 / den
            if n > 2:
                ref[np.arange(1, n - 1), np.arange(1, n - 1)] = (1 + rho**2) / den
            off = -rho / den
            ref += np.diag(np.repeat(off, n - 1), 1) + np.diag(np.repeat(off, n - 1), -1)
            max_regular = max(max_regular, float(np.max(np.abs(q - ref))))
    checks["regular_reduction_max_abs"] = max_regular

    thresholds = {
        "spd_inverse_identity_max_abs": 1e-11,
        "spd_solve_residual_max_abs": 1e-11,
        "iar1_precision_max_abs": 2e-8,
        "iar1_apply_max_abs": 2e-8,
        "iar1_profile_objective_max_abs": 2e-8,
        "regular_reduction_max_abs": 2e-11,
        "iar1_zero_precision_max_abs": 0.0,
        "iar1_zero_apply_max_abs": 0.0,
    }
    failures = {k: checks[k] for k, limit in thresholds.items() if float(checks[k]) > limit}
    checks["passed"] = not failures
    checks["failures"] = failures

    out = Path("validation/numerical_validation.json")
    out.write_text(json.dumps(checks, indent=2, sort_keys=True) + "\n")
    print(json.dumps(checks, indent=2, sort_keys=True))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
