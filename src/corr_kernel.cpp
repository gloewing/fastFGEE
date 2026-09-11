// Optional compiled Kronecker working-correlation inverse.
//
// This replaces .fgee_apply_kron_inverse() -- two reshapes, two aperm() copies
// and the two exact axis operators, each of which materialises several
// full-size temporaries -- with a single output allocation and two strided
// passes.  The R implementation remains the reference and the fallback; see
// .fgee_kron_kernel_try() in R/corr_kernel.R for the guard conditions and
// tests/testthat/test-corr-kernel.R for the agreement checks.
//
// Correlation codes: 0 = independent, 1 = exchangeable, 2 = AR(1) on a regular
// grid.  The cluster matrix is ordered with the functional index fastest, i.e.
// row f + n_fun * l, matching .fgee_tensor_grid_info()'s canonical ordering.
#include <Rcpp.h>
#include <cstring>
using namespace Rcpp;

// Exchangeable inverse along an axis of length n at stride s.  Sherman-Morrison:
// R^{-1} = I/(1-rho) - rho/((1-rho)(1+(n-1)rho)) J.  Each element needs only
// its own original value and the axis sum, so this is safe in place.
static inline void exch_axis(double* x, int n, int s, double rho) {
  if (n <= 1) return;
  const double d1 = 1.0 - rho, d2 = 1.0 + (n - 1) * rho;
  const double a = 1.0 / d1, b = rho / (d1 * d2);
  double sum = 0.0;
  for (int i = 0, k = 0; i < n; ++i, k += s) sum += x[k];
  const double sb = b * sum;
  for (int i = 0, k = 0; i < n; ++i, k += s) x[k] = a * x[k] - sb;
}

// Regular-grid AR(1) inverse along an axis: the tridiagonal precision matrix,
// applied in place with a two-element rolling buffer of the original values.
static inline void ar1_axis(double* x, int n, int s, double rho) {
  if (n <= 1) return;
  const double den = 1.0 - rho * rho, c = 1.0 + rho * rho;
  if (n == 2) {
    const double b0 = x[0], b1 = x[s];
    x[0] = (b0 - rho * b1) / den;
    x[s] = (b1 - rho * b0) / den;
    return;
  }
  double bprev = x[0];
  double bcur = x[s];
  x[0] = (bprev - rho * bcur) / den;
  for (int i = 1; i < n - 1; ++i) {
    const double bnext = x[(std::size_t)(i + 1) * s];
    x[(std::size_t) i * s] = (c * bcur - rho * (bprev + bnext)) / den;
    bprev = bcur;
    bcur = bnext;
  }
  x[(std::size_t)(n - 1) * s] = (bcur - rho * bprev) / den;
}

static inline void axis_apply(double* x, int n, int s, int code, double rho) {
  if (code == 1) exch_axis(x, n, s, rho);
  else if (code == 2) ar1_axis(x, n, s, rho);
}

// [[Rcpp::export(.fgee_kron_inverse_kernel)]]
NumericMatrix fgee_kron_inverse_kernel(const NumericMatrix& Q,
                                       int n_fun, int n_long,
                                       int code_fn, double rho_fn,
                                       int code_long, double rho_long) {
  const int m = Q.nrow(), q = Q.ncol();
  if (m != (double) n_fun * (double) n_long) {
    stop("Kronecker RHS has incompatible dimensions.");
  }
  if (code_fn < 0 || code_fn > 2 || code_long < 0 || code_long > 2) {
    stop("Unsupported correlation code.");
  }
  NumericMatrix O(m, q);
  if (m > 0 && q > 0) {
    std::memcpy(&O[0], &Q[0], (std::size_t) m * (std::size_t) q * sizeof(double));
  }
  double* o = (m > 0 && q > 0) ? &O[0] : 0;

  if (o && code_fn != 0 && n_fun > 1) {
    for (int c = 0; c < q; ++c) {
      double* col = o + (std::size_t) c * m;
      for (int l = 0; l < n_long; ++l) {
        axis_apply(col + (std::size_t) l * n_fun, n_fun, 1, code_fn, rho_fn);
      }
    }
  }
  if (o && code_long != 0 && n_long > 1) {
    for (int c = 0; c < q; ++c) {
      double* col = o + (std::size_t) c * m;
      for (int f = 0; f < n_fun; ++f) {
        axis_apply(col + f, n_long, n_fun, code_long, rho_long);
      }
    }
  }
  return O;
}
