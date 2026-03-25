#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================================
// Helper: cluster averaging
// CRITICAL: Assumes the data is pre-sorted by the 'grp' vector.
// ============================================================================
vec cluster_mean_cols_efficient(const mat& C, const ivec& grp) {
  int n = C.n_rows;
  int G = C.n_cols;
  vec out(G, fill::zeros);
  if (n == 0) return out;
  
  // Build cluster segments (start index + length) once
  std::vector<int> start;
  std::vector<int> len;
  start.reserve(n);
  len.reserve(n);
  
  int current = grp(0);
  int s = 0;
  int cnt = 1;
  for (int i = 1; i < n; ++i) {
    if (grp(i) == current) {
      cnt++;
    } else {
      start.push_back(s);
      len.push_back(cnt);
      s = i;
      cnt = 1;
      current = grp(i);
    }
  }
  start.push_back(s);
  len.push_back(cnt);
  
  int ncl = (int)start.size();
  
  // Column-major accumulation (cache-friendly)
  for (int g = 0; g < G; ++g) {
    const double* col = C.colptr(g);
    double grand = 0.0;
    
    for (int c = 0; c < ncl; ++c) {
      int st = start[c];
      int L  = len[c];
      
      double sum = 0.0;
      for (int i = 0; i < L; ++i) sum += col[st + i];
      
      grand += sum / (double)L;
    }
    
    out(g) = grand / (double)ncl;
  }
  
  return out;
}


// ============================================================================
// Link function: eta -> mu (inverse link)
// ============================================================================
mat apply_inverse_link(const mat& eta, const std::string& link, double clip_prob) {
  mat mu;
  
  if (link == "identity") {
    mu = eta;
  } else if (link == "log") {
    mu = exp(clamp(eta, -700.0, 700.0));
  } else if (link == "logit") {
    mu = 1.0 / (1.0 + exp(-clamp(eta, -700.0, 700.0)));
    if (clip_prob > 0) {
      mu.clamp(clip_prob, 1.0 - clip_prob);
    }
  } else if (link == "probit") {
    // optional clamp for numerical sanity (tails get clipped anyway)
    mat z = clamp(eta, -10.0, 10.0);
    
    // Φ(x) = 0.5 * erfc(-x / sqrt(2))
    mu = 0.5 * arma::erfc(-z / std::sqrt(2.0));
    
    if (clip_prob > 0) mu.clamp(clip_prob, 1.0 - clip_prob);
  } else if (link == "cloglog") {
    mat exp_eta = exp(clamp(eta, -700.0, 700.0));
    mu = 1.0 - exp(-exp_eta);
    if (clip_prob > 0) {
      mu.clamp(clip_prob, 1.0 - clip_prob);
    }
  } else if (link == "inverse") {
    mu = 1.0 / eta;
    mu.replace(datum::inf, 1e12); // Handle division by zero
    mu.replace(-datum::inf, -1e12);
  } else {
    stop("Unsupported link function: " + link);
  }
  
  return mu;
}

// ============================================================================
// Loss Functions (Robust, Vectorized, and Corrected)
// ============================================================================

vec gaussian_loss(const mat& mu, const vec& y, const ivec& grp, const std::string& loss_type) {
  // Corrected to be robust
  mat diff = mu.each_col() - y;
  mat C = square(diff);
  return cluster_mean_cols_efficient(C, grp);
}

vec binomial_loss(const mat& mu, const vec& y, const ivec& grp, const std::string& loss_type, double clip_prob) {
  mat C;
  if (loss_type == "brier") {
    // CORRECTED: Assign to the outer C, do not declare a new one.
    C = square(mu.each_col() - y);
  } else { // "nll"
    mat log_mu = log(mu);
    mat term1 = log_mu.each_col() % y;
    mat log_1_mu = log(1.0 - mu);
    mat term2 = log_1_mu.each_col() % (1.0 - y);
    C = -(term1 + term2);
  }
  return cluster_mean_cols_efficient(C, grp);
}

vec poisson_loss(const mat& mu, const vec& y, const ivec& grp, const std::string& loss_type) {
  mat C;
  if (loss_type == "brier") {
    // CORRECTED: Assign to the outer C, do not declare a new one.
    C = square(mu.each_col() - y);
  } else {  // "nll"
    mat mu_clipped = clamp(mu, 1e-10, 1e10);
    mat log_mu_term = log(mu_clipped);
    log_mu_term.each_col() %= y;
    C = mu_clipped - log_mu_term;
  }
  return cluster_mean_cols_efficient(C, grp);
}

vec gamma_loss(const mat& mu, const vec& y, const ivec& grp, const std::string& loss_type, double shape) {
  mat C;
  if (loss_type == "brier") {
    // CORRECTED: Assign to the outer C, do not declare a new one.
    C = square(mu.each_col() - y);
  } else { // "nll"
    mat mu_clipped = clamp(mu, 1e-10, 1e10);
    mat inv_mu = 1.0 / mu_clipped;
    inv_mu.each_col() %= y;
    C = log(mu_clipped) + inv_mu - 1.0;
  }
  return cluster_mean_cols_efficient(C, grp);
}


vec negbinom_loss(const mat& mu, const vec& y, const ivec& grp, double theta) {
  mat mu_clipped = clamp(mu, 1e-10, 1e10);
  mat mu_plus_theta = mu_clipped + theta;
  mat log_mu_plus_theta = log(mu_plus_theta);
  mat part1 = log_mu_plus_theta.each_col() % (y + theta);
  mat log_mu = log(mu_clipped);
  mat part2 = log_mu.each_col() % y;
  mat C = part1 - part2;
  return cluster_mean_cols_efficient(C, grp);
}

vec beta_loss(const mat& mu, const vec& y, const ivec& grp, double phi) {
  int G = mu.n_cols;
  mat C(mu.n_rows, G);
  vec y_clipped = clamp(y, 1e-6, 1.0 - 1e-6);
  vec log_y = log(y_clipped);
  vec log_1_y = log(1.0 - y_clipped);
  for (int g = 0; g < G; g++) {
    vec m = clamp(mu.col(g), 1e-6, 1.0 - 1e-6);
    vec a = m * phi;
    vec b = (1.0 - m) * phi;
    for (int i = 0; i < mu.n_rows; i++) {
      C(i, g) = R::lgammafn(a(i)) + R::lgammafn(b(i)) - R::lgammafn(a(i) + b(i))
      - (a(i) - 1.0) * log_y(i) 
      - (b(i) - 1.0) * log_1_y(i);
    }
  }
  return cluster_mean_cols_efficient(C, grp);
}

// ============================================================================
// Main dispatcher function
// ============================================================================
// [[Rcpp::export]]
vec fold_loss_from_eta_cluster_cpp(
    const mat& eta_mat,
    const vec& y,
    const ivec& grp,
    const std::string& family,
    const std::string& link,
    const std::string& loss,
    double clip_prob = 1e-6,
    double dispersion = 1.0
) {
  if ((int)eta_mat.n_rows != (int)y.n_elem) {
    stop("eta_mat rows must match length(y)");
  }
  
  mat mu;
  // For gaussian identity, mu is eta and we use it directly in the loss function
  if (family == "gaussian" && link == "identity") {
    mu = eta_mat;
  } else {
    mu = apply_inverse_link(eta_mat, link, clip_prob);
  }
  
  if (family == "gaussian") {
    return gaussian_loss(mu, y, grp, loss);
  } else if (family == "binomial" || family == "quasibinomial") {
    return binomial_loss(mu, y, grp, loss, clip_prob);
  } else if (family == "poisson" || family == "quasipoisson") {
    return poisson_loss(mu, y, grp, loss);
  } else if (family == "gamma") {
    return gamma_loss(mu, y, grp, loss, dispersion);
  } else if (family == "negbinomial" || family == "negative.binomial") {
    return negbinom_loss(mu, y, grp, dispersion);
  } else if (family == "beta") {
    return beta_loss(mu, y, grp, dispersion);
  } else {
    stop("Unsupported family: " + family);
  }
}

// ============================================================================
// Standalone wrapper for single fold testing
// ============================================================================
// [[Rcpp::export]]
NumericVector compute_fold_loss_cpp(
    NumericMatrix eta_mat_,
    NumericVector y_,
    IntegerVector grp_,
    std::string family,
    std::string link,
    std::string loss,
    double clip_prob = 1e-6,
    double dispersion = 1.0
) {
  mat eta_mat = as<mat>(eta_mat_);
  vec y = as<vec>(y_);
  ivec grp = as<ivec>(grp_);
  
  vec result = fold_loss_from_eta_cluster_cpp(
    eta_mat, y, grp, family, link, loss, clip_prob, dispersion
  );
  
  return wrap(result);
}

// ============================================================================
// Fast K-fold wrapper for all folds
// ============================================================================
// [[Rcpp::export]]
NumericMatrix fastkfold_loss_all_folds_cpp(
    const List& XX_list,
    const List& YY_list,
    const List& GRP_list,
    const arma::mat& beta,           // <- was NumericMatrix beta_
    const arma::cube& delta_arr,     // <- was arma::cube delta_arr (by value)
    const std::string& family,
    const std::string& link,
    const std::string& loss,
    double clip_prob = 1e-6,
    double dispersion = 1.0
) {
  int K = XX_list.size();
  int G = delta_arr.n_slices;
  int p = delta_arr.n_rows;
  
  arma::mat mse_mat(G, K);
  
  for (int kk = 0; kk < K; kk++) {
    
    // ---- X view (NO COPY) ----
    NumericMatrix Xr = XX_list[kk];
    arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
    
    // ---- y view (NO COPY) ----
    NumericVector yr = YY_list[kk];
    arma::vec y(yr.begin(), yr.size(), false);
    
    // ---- grp view (NO COPY) ----
    IntegerVector gr = GRP_list[kk];
    arma::ivec grp(gr.begin(), gr.size(), false);
    
    arma::mat Bhat(p, G);
    Bhat.each_col() = beta.col(kk);
    for (int g = 0; g < G; g++) {
      Bhat.col(g) += delta_arr.slice(g).col(kk);
    }
    
    arma::mat eta(X.n_rows, G);
    // arma::noalias(eta) = X * Bhat;
    eta = X * Bhat;
    
    arma::vec losses = fold_loss_from_eta_cluster_cpp(
      eta, y, grp, family, link, loss, clip_prob, dispersion
    );
    
    mse_mat.col(kk) = losses;
  }
  
  return wrap(mse_mat);
}

