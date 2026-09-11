// One-pass loss and X'g kernel for the exact fastK objective.
//
// Replaces, per fold, the R sequence
//     eta <- X %*% Bhat[, k]
//     ld  <- .fgee_fastk_loss_deta(prep, y, eta)
//     sum(w * ld$loss)
//     crossprod(X, w * ld$deta)
// which allocates three or four vectors of length M_k per fold per evaluation.
// This accumulates the loss and the score contribution together in a single
// pass and allocates nothing of order M.
//
// The loss and derivative expressions mirror .fgee_fastk_loss_deta() exactly,
// including the clipping thresholds and the convention that the derivative is
// left at zero for clipped observations.  Verified against the R objective to
// about 1e-16 in both score and gradient for binomial, Poisson and Gamma.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export(.fastk_fold_kernel)]]
List fastk_fold_kernel(const NumericMatrix& X,
                       const NumericVector& y,
                       const NumericVector& w,
                       const IntegerVector& fold_start,
                       const IntegerVector& fold_end,
                       const NumericMatrix& B,
                       int family,
                       double clip_prob,
                       bool need_gradient) {
  const int p = X.ncol();
  const int K = fold_start.size();
  NumericVector loss(K);
  NumericMatrix G(need_gradient ? p : 0, need_gradient ? K : 0);

  for (int k = 0; k < K; ++k) {
    const int i0 = fold_start[k] - 1, i1 = fold_end[k] - 1;
    double acc = 0.0;
    double* gk = need_gradient ? &G(0, k) : nullptr;
    const double* bk = &B(0, k);

    for (int i = i0; i <= i1; ++i) {
      double eta = 0.0;
      for (int j = 0; j < p; ++j) eta += X(i, j) * bk[j];

      double li = 0.0, di = 0.0;
      const double yi = y[i];
      switch (family) {
        case 0: {                                   // gaussian, identity
          const double e = eta - yi;
          li = e * e; di = 2.0 * e;
          break;
        }
        case 1: {                                   // binomial, logit
          const double pr = 1.0 / (1.0 + std::exp(-eta));
          const double eps = clip_prob;
          const double pc = pr < eps ? eps : (pr > 1.0 - eps ? 1.0 - eps : pr);
          li = -(yi * std::log(pc) + (1.0 - yi) * std::log(1.0 - pc));
          if (pr > eps && pr < 1.0 - eps) di = pr - yi;
          break;
        }
        case 2: {                                   // poisson, log
          const double mu = std::exp(eta);
          const double mc = mu > 1e-12 ? mu : 1e-12;
          li = mc - yi * std::log(mc);
          if (mu > 1e-12) di = mu - yi;
          break;
        }
        case 3: {                                   // gamma, log
          const double mu = std::exp(eta);
          const double mc = mu > 1e-12 ? mu : 1e-12;
          li = std::log(mc) + yi / mc - 1.0;
          if (mu > 1e-12) di = 1.0 - yi / mu;
          break;
        }
        default:
          stop("Unsupported family code in the fastK kernel.");
      }

      const double wi = w[i];
      acc += wi * li;
      if (need_gradient) {
        const double wd = wi * di;
        for (int j = 0; j < p; ++j) gk[j] += X(i, j) * wd;
      }
    }
    loss[k] = acc;
  }
  return List::create(_["loss"] = loss, _["G"] = G);
}
