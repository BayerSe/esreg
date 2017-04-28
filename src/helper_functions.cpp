#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List loop_esreg_covariance(arma::mat x, arma::colvec xq, arma::colvec xe,
                                 arma::colvec G1_prime_xq, arma::colvec G2_xe, arma::colvec G2_prime_xe,
                                 arma::colvec density, arma::colvec conditional_variance, double alpha) {
  int n = x.n_rows;
  int k = x.n_cols;

  arma::mat lambda_11 = arma::zeros<arma::mat>(k, k);
  arma::mat lambda_22 = arma::zeros<arma::mat>(k, k);
  arma::mat C_11 = arma::zeros<arma::mat>(k, k);
  arma::mat C_22 = arma::zeros<arma::mat>(k, k);
  arma::mat C_12 = arma::zeros<arma::mat>(k, k);

  arma::mat lambda = arma::zeros<arma::mat>(2*k, 2*k);
  arma::mat C = arma::zeros<arma::mat>(2*k, 2*k);
  arma::mat lambda_qr = arma::zeros<arma::mat>(k, k);
  arma::mat C_qr = arma::zeros<arma::mat>(k, k);

  arma::mat xi, xx;

  for (int i = 0; i < n; i++) {
    xi = x.row(i);
    xx = xi.t() * xi;
    lambda_11 += xx * (alpha * G1_prime_xq(i) + G2_xe(i)) / alpha * density(i);
    lambda_22 += xx * G2_prime_xe(i);
    C_11 += (1-alpha) / alpha * xx * pow((alpha * G1_prime_xq(i) + G2_xe(i)), 2);
    C_12 += (1-alpha) / alpha * xx * (alpha * G1_prime_xq(i) + G2_xe(i)) * G2_prime_xe(i) * (xq(i) - xe(i));
    C_22 += xx * pow(G2_prime_xe(i), 2) * (conditional_variance(i) / alpha + (1-alpha) / alpha * pow(xq(i) - xe(i), 2));

    lambda_qr += xx * (alpha * 1) / alpha * density(i);
    C_qr += (1-alpha) / alpha * xx * pow(alpha, 2);
  }

  C.submat(0, 0, k-1, k-1) = C_11;
  C.submat(0, k, k-1, 2*k-1) = C_12;
  C.submat(k, 0, 2*k-1, k-1) = C_12;
  C.submat(k, k, 2*k-1, 2*k-1) = C_22;

  lambda.submat(0, 0, k-1, k-1) = lambda_11;
  lambda.submat(k, k, 2*k-1, 2*k-1) = lambda_22;

  return Rcpp::List::create(Rcpp::Named("C") = C,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("C_qr") = C_qr,
                            Rcpp::Named("lambda_qr") = lambda_qr
  );
}


// [[Rcpp::export]]
arma::mat .esreg_twostep_covariance(arma::mat x, arma::colvec xq, arma::colvec xe,
                                    arma::colvec density, arma::colvec conditional_variance, double alpha) {
  try {
    int n = x.n_rows;
    int k = x.n_cols;

    // Define some 0-matrices
    arma::mat C_11 = arma::zeros<arma::mat>(k, k);
    arma::mat C_22 = arma::zeros<arma::mat>(k, k);

    arma::mat G_11 = arma::zeros<arma::mat>(k, k);
    arma::mat G_12 = arma::zeros<arma::mat>(k, k);
    arma::mat G_22 = arma::zeros<arma::mat>(k, k);

    arma::mat C = arma::zeros<arma::mat>(2*k, 2*k);
    arma::mat B = arma::zeros<arma::mat>(2*k, 2*k);

    arma::mat xi, xx;

    // Compute the matrix elements
    for (int i = 0; i < n; i++) {
      xi = x.row(i);
      xx = xi.t() * xi;

      C_11 += alpha * xx * conditional_variance(i);
      C_22 += alpha * (1 - alpha) * xx;

      G_11 += -alpha * xx;
      G_12 += xx * density(i) * (xe(i) - xq(i));
      G_22 += xx * density(i);
    }

    // Compute the averages
    C_11 /= n;
    C_22 /= n;
    G_11 /= n;
    G_12 /= n;
    G_22 /= n;

    // Fill the big matrices
    C.submat(0, 0, k-1, k-1) = C_11;
    C.submat(k, k, 2*k-1, 2*k-1) = C_22;

    B.submat(0, 0, k-1, k-1) = -inv(G_11);
    B.submat(0, k, k-1, 2*k-1) = -inv(G_11) * G_12 * inv(G_22);
    B.submat(k, k, 2*k-1, 2*k-1) = -inv(G_22);

    // Compute the covariance
    arma::mat cov = B * C * B.t();

    return cov / n;
  } catch(...) {
    Rcpp::stop("Cannot compute the two-step covariance!");
  }
}


// [[Rcpp::export]]
Rcpp::NumericMatrix stationary_bootstrap_indices(int n, double avg_block_size, int B) {
  Rcpp::NumericMatrix indices(n, B);
  int start_index, block_size;

  for (int b; b < B; b++) {
    int idx = 0;
    while (idx < n) {
      // Draw a random starting index
      start_index = Rcpp::as<int>(Rcpp::runif(1, 0, n));

      // Draw the block length
      block_size = 1 + Rcpp::as<int>(Rcpp::rgeom(1, 1 / avg_block_size));

      // Add the indices
      for (int i = start_index; i < start_index + block_size; i++) {
        // Wrap in a circle
        if (idx < n) {
          if (i < n) {
            indices(idx, b) = i;
          } else {
            indices(idx, b) = i - n;
          }
        }
        idx += 1;
      }
    }
  }
  return indices;
}



