#include <Rcpp.h>
using namespace Rcpp;

// One axis of sufficient statistics for the rank-one Gram fast path.
//
// Layout matches fgee_kron_inverse_kernel(): row index = f + n_fun * l, i.e.
// the functional index varies fastest within a cluster.
//
//   by_fn = true   ->  n_fun  x k, each functional point summed over the
//                      longitudinal index
//   by_fn = false  ->  n_long x k, each longitudinal block summed over the
//                      functional index
//
// Only the axis the caller needs is formed; see R/corr_gram.R.
//
// [[Rcpp::export(.fgee_gram_axis_sums)]]
NumericMatrix fgee_gram_axis_sums(const NumericMatrix& Q, int n_fun, int n_long,
                                  bool by_fn) {
  const R_xlen_t m = Q.nrow(), k = Q.ncol();
  if (n_fun <= 0 || n_long <= 0) stop("n_fun and n_long must be positive.");
  if (m != (R_xlen_t) n_fun * (R_xlen_t) n_long) {
    stop("Gram kernel RHS has incompatible dimensions.");
  }
  NumericMatrix S(by_fn ? n_fun : n_long, k);
  const double* q = &Q[0];
  for (R_xlen_t c = 0; c < k; ++c) {
    const double* col = q + c * m;
    if (by_fn) {
      double* s = &S[0] + c * (R_xlen_t) n_fun;
      for (int l = 0; l < n_long; ++l) {
        const double* b = col + (R_xlen_t) l * n_fun;
        for (int f = 0; f < n_fun; ++f) s[f] += b[f];
      }
    } else {
      double* s = &S[0] + c * (R_xlen_t) n_long;
      for (int l = 0; l < n_long; ++l) {
        const double* b = col + (R_xlen_t) l * n_fun;
        double acc = 0.0;
        for (int f = 0; f < n_fun; ++f) acc += b[f];
        s[l] = acc;
      }
    }
  }
  return S;
}
