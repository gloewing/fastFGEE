// Internal numerical kernels for fastFGEE.
//
// The irregular-time AR(1) implementation is an independent implementation of
// the Gaussian Markov precision and profiled likelihood described by
// Allevius (2018), "On the precision matrix of an irregularly sampled AR(1)
// process". No source code from the archived irregulAR1 package is used.
//
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>

#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

inline NumericMatrix checked_symmetric_copy(const NumericMatrix& x) {
  const int n = x.nrow();
  if (n < 1 || x.ncol() != n) {
    stop("x must be a non-empty square matrix");
  }

  NumericMatrix out(clone(x));
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if (!R_FINITE(out(i, j))) stop("x must contain only finite values");
    }
  }

  // Work with the symmetric part. The working breads are theoretically
  // symmetric, but floating-point accumulation can leave tiny asymmetry.
  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) {
      const double value = 0.5 * (out(i, j) + out(j, i));
      out(i, j) = value;
      out(j, i) = value;
    }
  }
  return out;
}

inline void chol_factor_upper(NumericMatrix& factor) {
  const La_INT n = static_cast<La_INT>(factor.nrow());
  La_INT info = 0;
  const char uplo = 'U';
  F77_CALL(dpotrf)(&uplo, &n, factor.begin(), &n, &info FCONE);
  if (info < 0) stop("LAPACK dpotrf received an invalid argument");
  if (info > 0) stop("matrix is not numerically positive definite");
}

inline void validate_rho(const double rho) {
  if (!R_FINITE(rho) || rho < 0.0 || rho >= 1.0) {
    stop("rho must be finite and satisfy 0 <= rho < 1");
  }
}

inline void validate_times(const NumericVector& time, const bool require_two) {
  const R_xlen_t n = time.size();
  if (n < (require_two ? 2 : 1)) {
    stop(require_two ? "at least two observation times are required" :
                       "at least one observation time is required");
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!R_FINITE(time[i])) stop("observation times must be finite");
    if (i > 0 && !(time[i] > time[i - 1])) {
      stop("observation times must be strictly increasing");
    }
  }
}

struct Transition {
  double a;
  double den;
};

inline Transition transition(const double rho, const double gap) {
  if (!R_FINITE(gap) || gap <= 0.0) {
    stop("irregular AR(1) time gaps must be positive and finite");
  }
  // rho = 0 is the independence boundary and is valid for every positive gap.
  if (rho == 0.0) {
    Transition out = {0.0, 1.0};
    return out;
  }
  const double log_a = std::log(rho) * gap;
  const double a = std::exp(log_a); // underflow to zero is a valid limit
  const double den = -std::expm1(2.0 * log_a); // 1 - a^2, stably
  if (!R_FINITE(a) || a < 0.0 || a >= 1.0 ||
      !R_FINITE(den) || den <= 0.0) {
    stop("invalid irregular AR(1) transition; check rho and time gaps");
  }
  Transition out = {a, den};
  return out;
}

inline void precision_bands(const NumericVector& time,
                            const double rho,
                            NumericVector& diagonal,
                            NumericVector& off_diagonal) {
  validate_rho(rho);
  validate_times(time, false);
  const R_xlen_t n = time.size();
  diagonal = NumericVector(n, 0.0);
  off_diagonal = NumericVector(n > 1 ? n - 1 : 0, 0.0);

  if (n == 1) {
    diagonal[0] = 1.0;
    return;
  }

  std::vector<double> a(static_cast<std::size_t>(n - 1), 0.0);
  std::vector<double> den(static_cast<std::size_t>(n - 1), 0.0);
  for (R_xlen_t j = 1; j < n; ++j) {
    const Transition tr = transition(rho, time[j] - time[j - 1]);
    a[static_cast<std::size_t>(j - 1)] = tr.a;
    den[static_cast<std::size_t>(j - 1)] = tr.den;
  }

  diagonal[0] = 1.0 / den[0];
  for (R_xlen_t j = 1; j < n - 1; ++j) {
    const double left = 1.0 / den[static_cast<std::size_t>(j - 1)];
    const double right_a = a[static_cast<std::size_t>(j)];
    const double right = right_a * right_a /
      den[static_cast<std::size_t>(j)];
    diagonal[j] = left + right;
  }
  diagonal[n - 1] = 1.0 / den[static_cast<std::size_t>(n - 2)];

  for (R_xlen_t j = 0; j < n - 1; ++j) {
    off_diagonal[j] = -a[static_cast<std::size_t>(j)] /
      den[static_cast<std::size_t>(j)];
  }
}

} // namespace

// Invert a real symmetric positive-definite matrix with LAPACK dpotrf/dpotri.
// [[Rcpp::export]]
NumericMatrix fgee_sympd_inverse_cpp(const NumericMatrix& x) {
  NumericMatrix ans = checked_symmetric_copy(x);
  const La_INT n = static_cast<La_INT>(ans.nrow());
  chol_factor_upper(ans);

  La_INT info = 0;
  const char uplo = 'U';
  F77_CALL(dpotri)(&uplo, &n, ans.begin(), &n, &info FCONE);
  if (info < 0) stop("LAPACK dpotri received an invalid argument");
  if (info > 0) stop("the Cholesky factor is numerically singular");

  for (int j = 0; j < n; ++j) {
    for (int i = j + 1; i < n; ++i) ans(i, j) = ans(j, i);
  }
  return ans;
}

// Solve A X = B for real symmetric positive-definite A with dpotrf/dpotrs.
// [[Rcpp::export]]
NumericMatrix fgee_sympd_solve_cpp(const NumericMatrix& x,
                                  const NumericMatrix& b) {
  NumericMatrix factor = checked_symmetric_copy(x);
  const La_INT n = static_cast<La_INT>(factor.nrow());
  if (b.nrow() != n) stop("nrow(b) must equal nrow(x)");

  NumericMatrix ans(clone(b));
  const La_INT nrhs = static_cast<La_INT>(ans.ncol());
  for (int j = 0; j < ans.ncol(); ++j) {
    for (int i = 0; i < ans.nrow(); ++i) {
      if (!R_FINITE(ans(i, j))) stop("b must contain only finite values");
    }
  }

  chol_factor_upper(factor);
  La_INT info = 0;
  const char uplo = 'U';
  F77_CALL(dpotrs)(&uplo, &n, &nrhs, factor.begin(), &n,
                   ans.begin(), &n, &info FCONE);
  if (info < 0) stop("LAPACK dpotrs received an invalid argument");
  return ans;
}

// Return the diagonal and first off-diagonal of the exact Markov precision for
// Corr(Z_j, Z_k) = rho^abs(t_j - t_k), including rho = 0 as the
// working-independence boundary.
// [[Rcpp::export]]
List fgee_iar1_precision_bands_cpp(const NumericVector& time,
                                  const double rho) {
  NumericVector diagonal;
  NumericVector off_diagonal;
  precision_bands(time, rho, diagonal, off_diagonal);
  return List::create(_["diagonal"] = diagonal,
                      _["off_diagonal"] = off_diagonal);
}

// Dense reference form of the irregular AR(1) precision.
// [[Rcpp::export]]
NumericMatrix fgee_iar1_precision_cpp(const NumericVector& time,
                                     const double rho) {
  NumericVector diagonal;
  NumericVector off_diagonal;
  precision_bands(time, rho, diagonal, off_diagonal);
  const R_xlen_t n = time.size();
  NumericMatrix out(n, n);
  for (R_xlen_t j = 0; j < n; ++j) out(j, j) = diagonal[j];
  for (R_xlen_t j = 0; j + 1 < n; ++j) {
    out(j, j + 1) = off_diagonal[j];
    out(j + 1, j) = off_diagonal[j];
  }
  return out;
}

// Apply the exact tridiagonal precision to one or more right-hand sides.
// [[Rcpp::export]]
NumericMatrix fgee_iar1_apply_precision_cpp(const NumericMatrix& rhs,
                                           const NumericVector& time,
                                           const double rho) {
  const R_xlen_t n = time.size();
  if (rhs.nrow() != n) stop("nrow(rhs) must equal length(time)");

  NumericVector diagonal;
  NumericVector off_diagonal;
  precision_bands(time, rho, diagonal, off_diagonal);
  const int q = rhs.ncol();
  NumericMatrix out(n, q);
  for (int k = 0; k < q; ++k) {
    for (R_xlen_t j = 0; j < n; ++j) {
      const double current = rhs(j, k);
      if (!R_FINITE(current)) stop("rhs must contain only finite values");
      double value = diagonal[j] * current;
      if (j > 0) value += off_diagonal[j - 1] * rhs(j - 1, k);
      if (j + 1 < n) value += off_diagonal[j] * rhs(j + 1, k);
      out(j, k) = value;
    }
  }
  return out;
}

// Profiled Gaussian negative twice log likelihood, up to an additive constant,
// for a zero-mean irregularly sampled AR(1) series. The innovation variance is
// profiled out; multiplication by a positive constant does not change rho.
// [[Rcpp::export]]
double fgee_iar1_profile_nll_cpp(const NumericVector& residual,
                                const NumericVector& time,
                                const double rho) {
  validate_rho(rho);
  validate_times(time, true);
  const R_xlen_t n = residual.size();
  if (n != time.size()) stop("residual and time must have equal length");
  for (R_xlen_t j = 0; j < n; ++j) {
    if (!R_FINITE(residual[j])) stop("residuals must be finite");
  }

  double quadratic = residual[0] * residual[0];
  double log_determinant = 0.0;
  for (R_xlen_t j = 1; j < n; ++j) {
    const Transition tr = transition(rho, time[j] - time[j - 1]);
    const double innovation = residual[j] - tr.a * residual[j - 1];
    quadratic += innovation * innovation / tr.den;
    log_determinant += std::log(tr.den);
  }
  if (!(quadratic > 0.0) || !R_FINITE(quadratic)) {
    return std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(n) *
    std::log(quadratic / static_cast<double>(n)) + log_determinant;
}
