#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @keywords internal
// [[Rcpp::export]]
arma::mat lambda_matrix_loop(
    arma::mat xq, arma::mat xe, arma::vec xbq, arma::vec xbe,
    arma::vec G1_prime_xq, arma::vec G2_xe, arma::vec G2_prime_xe, arma::vec G2_prime_prime_xe,
    arma::vec density, bool include_misspecification_terms, double alpha) {
  int n = xq.n_rows;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Define some 0-matrices
  arma::mat lambda_11 = arma::zeros<arma::mat>(kq, kq);
  arma::mat lambda_12 = arma::zeros<arma::mat>(kq, ke);
  arma::mat lambda_22 = arma::zeros<arma::mat>(ke, ke);
  arma::mat lambda = arma::zeros<arma::mat>(kq+ke, kq+ke);
  arma::mat xqi, xei, xxq, xxe, xxqe;
  double yi, xbqi, xbei, hit;

  // Compute the matrix elements
  for (int i = 0; i < n; i++) {
    xqi = xq.row(i);
    xei = xe.row(i);
    xxq = xqi.t() * xqi;
    xxe = xei.t() * xei;
    xxqe = xqi.t() * xei;
    xbqi = xbq(i);
    xbei = xbe(i);
    hit = yi <= xbqi;

    if (include_misspecification_terms) {
      lambda_11 += xxq * (G1_prime_xq(i) + G2_xe(i) / alpha) * density(i);
      lambda_12 += xxqe * G2_prime_xe(i) * (hit - alpha) / alpha;
      lambda_22 += xxe * G2_prime_xe(i) + xxe * G2_prime_prime_xe(i) *(xbei - yi * hit / alpha + xbqi * (hit - alpha) / alpha);
    } else {
      lambda_11 += xxq * (G1_prime_xq(i) + G2_xe(i) / alpha) * density(i);
      lambda_22 += xxe * G2_prime_xe(i);
    }
  }

  // Fill the matrices
  lambda.submat(0, 0, kq-1, kq-1) = lambda_11 / n;
  lambda.submat(0, kq, kq-1, kq+ke-1) = lambda_12 / n;
  lambda.submat(kq, 0, kq+ke-1, kq-1) = lambda_12.t() / n;
  lambda.submat(kq, kq, kq+ke-1, kq+ke-1) = lambda_22 / n;

  return lambda;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat sigma_matrix_calculcated(
    arma::mat xq, arma::mat xe, arma::colvec xbq, arma::colvec xbe,
    arma::colvec G1_prime_xq, arma::colvec G2_xe, arma::colvec G2_prime_xe,
    arma::colvec conditional_variance, double alpha) {
  int n = xq.n_rows;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Define some 0-matrices
  arma::mat sigma_11 = arma::zeros<arma::mat>(kq, kq);
  arma::mat sigma_12 = arma::zeros<arma::mat>(ke, kq);
  arma::mat sigma_22 = arma::zeros<arma::mat>(ke, ke);
  arma::mat sigma = arma::zeros<arma::mat>(kq+ke, kq+ke);
  arma::mat xqi, xei, xxq, xxe, xxeq;

  // Compute the matrix elements
  for (int i = 0; i < n; i++) {
    xqi = xq.row(i);
    xei = xe.row(i);
    xxq = xqi.t() * xqi;
    xxe = xei.t() * xei;
    xxeq = xei.t() * xqi;

    sigma_11 += (1-alpha)/alpha * xxq * pow(alpha*G1_prime_xq(i) + G2_xe(i), 2);
    sigma_12 += (1-alpha)/alpha * xxeq * (xbq(i) - xbe(i)) *
      (alpha*G1_prime_xq(i) + G2_xe(i)) * G2_prime_xe(i);
    sigma_22 += xxe * pow(G2_prime_xe(i), 2) * (conditional_variance(i)/alpha +
      (1-alpha)/alpha * pow(xbq(i) - xbe(i), 2));
  }

  // Fill the matrices
  sigma.submat(0, 0, kq-1, kq-1) = sigma_11 / n;
  sigma.submat(kq, 0, kq+ke-1, kq-1) = sigma_12 / n;
  sigma.submat(0, kq, kq-1, kq+ke-1) = sigma_12.t() / n;
  sigma.submat(kq, kq, kq+ke-1, kq+ke-1) = sigma_22 / n;

  return sigma;
}



//' @keywords internal
// [[Rcpp::export]]
arma::mat estimating_function_loop(
    arma::vec y, arma::mat xq, arma::mat xe, arma::colvec xbq, arma::colvec xbe,
    arma::colvec G1_prime_xq, arma::colvec G2_xe, arma::colvec G2_prime_xe, double alpha) {
  int n = xq.n_rows;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Initialize variables
  arma::mat psi = arma::zeros<arma::mat>(n, kq + ke);
  double yi, xbqi, xbei, hit;
  arma::mat xqi, xei, xxq, xxe;

  // Compute the matrix elements
  for (int i = 0; i < n; i++) {
    yi = y(i);
    xqi = xq.row(i);
    xei = xe.row(i);
    xxq = xqi.t() * xqi;
    xxe = xei.t() * xei;
    xbqi = xbq(i);
    xbei = xbe(i);

    // Hit variable
    hit = yi <= xbqi;

    // Fill the matrix
    psi.submat(i, 0, i, kq-1) = xqi * (G1_prime_xq(i) + G2_xe(i)/alpha) * (hit - alpha);
    psi.submat(i, kq, i, kq+ke-1) = xei *  G2_prime_xe(i) * (xbei - xbqi + (xbqi - yi) * hit / alpha);
  }

  return psi;
}
