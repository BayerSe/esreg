#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Specification Function
//' @description G1
//' @param z Data
//' @param type Choice of the G1 function:
//' \itemize{
//'   \item 1: G1(z) = z
//'   \item 2: G1(z) = 0
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G1_fun(double z, int type) {
  double out;
  if (type == 1) {
    out = z;
  } else if (type == 2) {
    out = 0;
  } else {
    Rcpp::stop("type not in 1, 2!");
  }
  return out;
}

//' @title Specification Function
//' @description G1_prime
//' @param z Data
//' @param type Choice of the G1_prime function:
//' \itemize{
//'   \item 1: G1_prime(z) = 1
//'   \item 2: G1_prime(z) = 0
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G1_prime_fun(double z, int type) {
  double out;
  if (type == 1) {
    out = 1;
  } else if (type == 2) {
    out = 0;
  } else {
    Rcpp::stop("type not in 1, 2!");
  }
  return out;
}


//' @title Specification Function
//' @description G2_curly
//' @param z Data
//' @param type Choice of the G2_curly function:
//' \itemize{
//'   \item 1: -log(-z), z < 0
//'   \item 2: -sqrt(-z), z < 0
//'   \item 3: -1/z, z < 0
//'   \item 4: log(1 + exp(z))
//'   \item 5: exp(z)
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G2_curly_fun(double z, int type) {
  double out;
  if ((type == 1) | (type == 2) | (type == 3)) {
    if (z >= 0) {
      Rcpp::warning("z can not be positive for type 1, 2, 3!");
      return NA_REAL;
    }
  }
  if (type == 1) {
    out = -log(-z);
  } else if (type == 2) {
    out = -sqrt(-z);
  } else if (type == 3) {
    out = -1/z;
  } else if (type == 4) {
    out = log1p(exp(z));
  } else if (type == 5) {
    out = exp(z);
  } else {
    Rcpp::stop("type not in 1, 2, 3, 4, 5!");
  }
  return out;
}


//' @title Specification Function
//' @description G2
//' @param z Data
//' @param type Choice of the G2 function:
//' \itemize{
//'   \item 1: -1/z, z < 0
//'   \item 2: 0.5/sqrt(-z), z < 0
//'   \item 3: 1/z^2, z < 0
//'   \item 4: 1 / (1 + exp(-z))
//'   \item 5: exp(z)
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G2_fun(double z, int type) {
  double out;
  if ((type == 1) | (type == 2) | (type == 3)) {
    if (z >= 0) {
      Rcpp::warning("z can not be positive for type 1, 2, 3!");
      return NA_REAL;
    }
  }
  if (type == 1) {
    out = -1/z;
  } else if (type == 2) {
    out = 0.5/sqrt(-z);
  } else if (type == 3) {
    out = 1/pow(z, 2);
  } else if (type == 4) {
    out = 1/(1 + exp(-z));
  } else if (type == 5) {
    out = exp(z);
  } else {
    Rcpp::stop("type not in 1, 2, 3, 4, 5!");
  }
  return out;
}


//' @title Specification Function
//' @description G2_prime
//' @param z Data
//' @param type Choice of the G2_prime function:
//' \itemize{
//'   \item 1: 1/z^2, z < 0
//'   \item 2: 0.25 / (-z)^(3/2), z < 0
//'   \item 3: -2/z^3, z < 0
//'   \item 4: exp(z) / (1 + exp(z))^2
//'   \item 5: exp(z)
//' }
//' @keywords internal
//' @export
// [[Rcpp::export]]
double G2_prime_fun(double z, int type) {
  double out;
  if ((type == 1) | (type == 2) | (type == 3)) {
    if (z >= 0) {
      Rcpp::warning("z can not be positive for type 1, 2, 3!");
      return NA_REAL;
    }
  }
  if (type == 1) {
    out = 1/pow(z, 2);
  } else if (type == 2) {
    out = 0.25/pow(-z, 1.5);
  } else if (type == 3) {
    out = -2/pow(z, 3);
  } else if (type == 4) {
    out =  exp(z) / pow(1 + exp(z), 2);
  } else if (type == 5) {
    out = exp(z);
  } else {
    Rcpp::stop("type not in 1, 2, 3, 4, 5!");
  }
  return out;
}

//' @title Vectorized call to the G1 / G2 functions
//' @description Vectorized call to the G1 / G2 functions
//' @param z Vector
//' @param g String, either G1, G1_prime, G2_curly, G2 or G2_curly
//' @param type Numeric, for G1: 1-2; G2: 1-5
//' (see \link{G1_fun}, \link{G1_prime_fun}, \link{G2_curly_fun}, \link{G2_fun}, \link{G2_prime_fun})
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector G_vec(Rcpp::NumericVector z, Rcpp::String g, int type) {
  int n = z.size();
  Rcpp::NumericVector out(n);
  if (g == "G1") {
    for (int i = 0; i < n; i++) out[i] = G1_fun(z[i], type);
  } else if (g == "G1_prime") {
    for (int i = 0; i < n; i++) out[i] = G1_prime_fun(z[i], type);
  } else if (g == "G2_curly") {
    for (int i = 0; i < n; i++) out[i] = G2_curly_fun(z[i], type);
  } else if (g == "G2") {
    for (int i = 0; i < n; i++) out[i] = G2_fun(z[i], type);
  } else if (g == "G2_prime") {
    for (int i = 0; i < n; i++) out[i] = G2_prime_fun(z[i], type);
  } else {
    Rcpp::stop("Non supported G-function!");
  }
  return out;
}


//' @title Joint (VaR, ES) loss for a linear predictor
//' @description Returns the loss for the parameter vector b
//' @param b Parameter vector
//' @param y Vector of dependent data
//' @param xq Matrix of covariates for the quantile part
//' @param xe Matrix of covariates for the expected shortfall part
//' @param alpha Probability level
//' @param g1 1, 2 (see \link{G1_fun})
//' @param g2 1, 2, 3, 4, 5 (see \link{G2_curly_fun}, \link{G2_fun})
//' @importFrom Rcpp sourceCpp
//' @useDynLib esreg
//' @keywords internal
//' @export
// [[Rcpp::export]]
double esr_rho_lp(const arma::colvec& b, const arma::colvec& y,
                  const arma::mat& xq, const arma::mat& xe,
                  double alpha, int g1=2L, int g2=1L) {

  int n = xq.n_rows;                          // Number of observations
  int kq = xq.n_cols;                         // Number of q-parameters
  int ke = xe.n_cols;                         // Number of es-parameters
  arma::colvec bq = b.subvec(0, kq-1);        // Quantile parameters
  arma::colvec be = b.subvec(kq, kq+ke-1);    // Expected shortfall parameters

  // Initialize variables
  double yi, xbq, xbe, h, out = 0;
  arma::mat xqi;
  arma::mat xei;

  // Compute the loss
  for (int i = 0; i < n; i++) {
    yi = y(i);
    xqi = xq.row(i).t();
    xei = xe.row(i).t();
    xbq = as_scalar(xqi.t() * bq);
    xbe = as_scalar(xei.t() * be);

    // Check the shortfall
    if (((g2 == 1) | (g2 == 2) | (g2 == 3)) & (xbe >= -0.01)) {
      Rcpp::warning("x'b_e can not be positive for g2 1, 2, 3!");
      return NA_REAL;
    }

    // Compute the loss
    h = yi <= xbq;  // Hit variable
    out += (h - alpha) * G1_fun(xbq, g1) - h * G1_fun(yi, g1) +
      G2_fun(xbe, g2) * (xbe - xbq + (xbq - yi) * h / alpha) -
      G2_curly_fun(xbe, g2);
  }

  return out / n;
}


//' @keywords internal
// [[Rcpp::export]]
arma::mat sigma_matrix(const Rcpp::List & object) {
  // Extract quantities from list input
  arma::vec bq = Rcpp::as<arma::vec>(object["coefficients_q"]);
  arma::vec be = Rcpp::as<arma::vec>(object["coefficients_e"]);

  arma::vec y = Rcpp::as<arma::vec>(object["y"]);
  arma::mat xq = Rcpp::as<arma::mat>(object["xq"]);
  arma::mat xe = Rcpp::as<arma::mat>(object["xe"]);

  double alpha = Rcpp::as<double>(object["alpha"]);
  int g1 = Rcpp::as<double>(object["g1"]);
  int g2 = Rcpp::as<double>(object["g2"]);

  // Extract dimensions
  int n = y.n_elem;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Initialize variables
  double yi, xbq, xbe, h;
  arma::mat xqi, xei;
  arma::mat psi = arma::zeros<arma::mat>(n, kq + ke);

  // Transform input variables
  if (((g2 == 1) | (g2 == 2) | (g2 == 3))) {
    double max_y = max(y);
    y -= max_y;
    bq(0) -= max_y;
    be(0) -= max_y;
  }

  // Loop over the observations
  for (int i = 0; i < n; i++) {
    yi = y(i);
    xqi = xq.row(i);
    xei = xe.row(i);
    xbq = as_scalar(xqi * bq);
    xbe = as_scalar(xei * be);

    // Check the shortfall
    if (((g2 == 1) | (g2 == 2) | (g2 == 3)) & (xbe >= 0)) {
      Rcpp::warning("x'b_e can not be positive for g2 1, 2, 3!");
      return arma::mat(n, kq+ke).fill(NA_REAL);
    }

    // Hit variable
    h = yi <= xbq;

    // Fill the matrix
    psi.submat(i, 0, i, kq-1) = xqi * (G1_prime_fun(xbq, g1) + G2_fun(xbe, g2)/alpha) * (h - alpha);
    psi.submat(i, kq, i, kq+ke-1) = xei *  G2_prime_fun(xbe, g2) * (xbe - xbq + (xbq - yi) * h / alpha);
  }

  // Compute the crossproduct
  arma::mat sigma = psi.t() * psi / n;

  return sigma;
}


//' @keywords internal
// [[Rcpp::export]]
arma::mat lambda_matrix(const Rcpp::List & object) {
  // Extract quantities from list input
  arma::vec bq = Rcpp::as<arma::vec>(object["coefficients_q"]);
  arma::vec be = Rcpp::as<arma::vec>(object["coefficients_e"]);

  arma::vec y = Rcpp::as<arma::vec>(object["y"]);
  arma::mat xq = Rcpp::as<arma::mat>(object["xq"]);
  arma::mat xe = Rcpp::as<arma::mat>(object["xe"]);

  double alpha = Rcpp::as<double>(object["alpha"]);
  int g1 = Rcpp::as<double>(object["g1"]);
  int g2 = Rcpp::as<double>(object["g2"]);

  // Extract dimensions
  int n = y.n_elem;
  int kq = xq.n_cols;
  int ke = xe.n_cols;

  // Define some 0-matrices
  arma::mat lambda_11 = arma::zeros<arma::mat>(kq, kq);
  arma::mat lambda_12 = arma::zeros<arma::mat>(kq, ke);
  arma::mat lambda_22 = arma::zeros<arma::mat>(ke, ke);

  arma::mat lambda = arma::zeros<arma::mat>(kq+ke, kq+ke);

  arma::mat xqi, xei, xxq, xxe, xxqe;
  double yi, xbqi, xbei, h, hit;
  double ct = pow(n, -1.0/3.0);
  //Rcpp::Rcout << "The value is " << ct << std::endl;

  // Transform input variables
  if (((g2 == 1) | (g2 == 2) | (g2 == 3))) {
    double max_y = max(y);
    y -= max_y;
    bq(0) -= max_y;
    be(0) -= max_y;
  }

  // Compute the matrix elements
  for (int i = 0; i < n; i++) {
    yi = y(i);

    xqi = xq.row(i);
    xei = xe.row(i);

    xxq = xqi.t() * xqi;
    xxe = xei.t() * xei;
    xxqe = xqi.t() * xei;

    xbqi = as_scalar(xqi * bq);
    xbei = as_scalar(xei * be);

    h = std::abs(yi - xbqi) <= ct;
    hit = yi <= xbqi;

    lambda_11 += -1/alpha * xxq / xbei * h / (2*ct);
    lambda_12 += xxqe / pow(xbei, 2) * (hit - alpha);
    lambda_22 += xxe / pow(xbei, 2);// - 2 * xxe / pow(xbei, 3) * (xbei - yi/alpha*hit + xbqi * (hit-alpha)/alpha);
  }

  // Fill the matrices
  lambda.submat(0, 0, kq-1, kq-1) = lambda_11 / n;
  // lambda.submat(0, kq, kq-1, kq+ke-1) = lambda_12 / n;
  // lambda.submat(kq, 0, kq+ke-1, kq-1) = lambda_12.t() / n;
  lambda.submat(kq, kq, kq+ke-1, kq+ke-1) = lambda_22 / n;

  return lambda;
}
