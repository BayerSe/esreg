#' Estimated asymptotic covariance for the joint estimator
#'
#' @param fit fit An object from calling esreg()
#' @param sparsity Sparsity estimator
#' \itemize{
#'    \item iid - Piecewise linear interpolation of the distribution
#'    \item nid - Hendricks and Koenker sandwich (two additional quantile regressions)
#' }
#' @param cond_var Conditional truncated variance estimator
#' \itemize{
#'    \item ind Variance over all negative residuals
#'    \item scl_N Scaling with the Normal distribution
#'    \item scl_t Scaling with the t-distribution
#' }
#' @param bandwidth_type Bofinger, Chamberlain or Hall-Sheather
#' @export
esreg_covariance <- function(fit, sparsity = "iid", cond_var = "ind",
                             bandwidth_type = "Hall-Sheather") {
  if (!methods::is(fit, "esreg"))
    stop("This is not a esreg object!")
  if (!(sparsity %in% c("iid", "nid")))
    stop("sparsity can be iid or nid")
  if (!(cond_var %in% c('ind', 'scl_N', 'scl_t')))
    stop('cond_var can be ind, scl_N or scl_t')
  if (!(bandwidth_type %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
    stop("bandwidth_type can be Bofinger, Chamberlain or Hall-Sheather")

  # Extract some elements from the esreg fit object
  y <- fit$y
  x <- fit$x
  alpha <- fit$alpha
  par_q <- fit$par_q
  par_e <- fit$par_e

  # Transform the data and coefficients
  if (fit$g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
    par_q[1] <- par_q[1] - max_y
    par_e[1] <- par_e[1] - max_y
  }

  # Precompute some quantities
  xq <- as.numeric(x %*% par_q)
  xe <- as.numeric(x %*% par_e)
  u <- as.numeric(y - xq)
  n <- nrow(x)
  k <- ncol(x)

  # Check the methods in case of sample quantile / es
  if ((k == 1) & sparsity != "iid") {
    warning("Changed sparsity estimation to iid!")
    sparsity <- "iid"
  }
  if ((k == 1) & cond_var != "ind") {
    warning("Changed condittional truncated variance estimation to nid!")
    cond_var <- "ind"
  }

  # Density quantile function
  dens <- density_quantile_function(y = y, x = x, u = u, alpha = alpha,
                                    sparsity = sparsity, bandwidth_type = bandwidth_type)

  # Truncated conditional variance
  cv <- conditional_truncated_variance(y = y, x = x, u = u, approach = cond_var)

  # Evaluate G1 / G2 functions
  G1_prime_xq <- G_vec(z = xq, g = 'G1_prime', type = fit$g1)
  G2_xe <- G_vec(z = xe, g = 'G2', type = fit$g2)
  G2_prime_xe <- G_vec(z = xe, g = 'G2_prime', type = fit$g2)

  # Compute and return the covariance matrix
  cov <- l_esreg_covariance(x = x, xq = xq, xe = xe, alpha = alpha,
                            G1_prime_xq = G1_prime_xq,
                            G2_xe = G2_xe, G2_prime_xe = G2_prime_xe,
                            density = dens, conditional_variance = cv)
  cov
}


#' Estimated asymptotic covariance for the two-step estimator
#'
#' @param fit fit An object from calling esreg_twostep()
#' @param sparsity Sparsity estimator
#' \itemize{
#'    \item iid - Piecewise linear interpolation of the distribution
#'    \item nid - Hendricks and Koenker sandwich (two additional quantile regressions)
#' }
#' @param cond_var Conditional truncated variance estimator
#' \itemize{
#'    \item ind Variance over all negative residuals
#'    \item scl_N Scaling with the Normal distribution
#'    \item scl_t Scaling with the t-distribution
#' }
#' @param bandwidth_type Bofinger, Chamberlain or Hall-Sheather
#' @export
esreg_twostep_covariance <- function(fit, sparsity = "iid", cond_var = "ind",
                                     bandwidth_type = "Hall-Sheather") {
  if (!methods::is(fit, "esreg_twostep"))
    stop("This is not a esreg_twostep object!")
  if (!(sparsity %in% c("iid", "nid")))
    stop("sparsity can be iid or nid")
  if (!(cond_var %in% c('ind', 'scl_N', 'scl_t')))
    stop('cond_var can be ind, scl_N or scl_t')
  if (!(bandwidth_type %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
    stop("bandwidth_type can be Bofinger, Chamberlain or Hall-Sheather")

  # Extract some elements from the esreg_twostep fit object
  y <- fit$y
  x <- fit$x
  alpha <- fit$alpha
  par_q <- fit$par_q
  par_e <- fit$par_e

  # Precompute some quantities
  xq <- as.numeric(x %*% par_q)
  xe <- as.numeric(x %*% par_e)
  u <- as.numeric(y - xq)
  n <- nrow(x)
  k <- ncol(x)

  # Check the methods in case of sample quantile / es
  if ((k == 1) & sparsity != "iid") {
    warning("Changed sparsity estimation to iid!")
    sparsity <- "iid"
  }
  if ((k == 1) & cond_var != "ind") {
    warning("Changed condittional truncated variance estimation to nid!")
    cond_var <- "ind"
  }

  # Density quantile function
  dens <- density_quantile_function(y = y, x = x, u = u, alpha = alpha,
                                    sparsity = sparsity, bandwidth_type = bandwidth_type)

  # Truncated conditional variance
  cv <- conditional_truncated_variance(y = y, x = x, u = u, approach = cond_var)

  # Compute and return the covariance matrix
  cov <- l_esreg_twostep_covariance(x = x, xq = xq, xe = xe, alpha = alpha,
                                    density = dens, conditional_variance = cv)

  # Flip blockwise such that it matches the rest of the esreg package
  cov <- rbind(cbind(cov[(k+1):(2*k), (k+1):(2*k)], cov[1:k, (k+1):(2*k)]),
               cbind(cov[(k+1):(2*k), 1:k], cov[1:k, 1:k]))
}


#' Bootstrapped Covariance
#'
#' Estimate the variance-covariance matrix via bootstrapping
#'
#' @param fit esreg object
#' @param B Number of bootstrap iterations
#' @param bootstrap_method iid or stationary
#' @param block_length Average block length for the stationary bootstrap
#' @references Politis & Romano (1994)
#' @export
esreg_covariance_boot <- function(fit, B = 1000, bootstrap_method = "iid", block_length = NULL) {

  # Draw the bootstrap indices
  n <- length(fit$y)
  if (bootstrap_method == "iid") {
    idx <- matrix(sample(1:n, size = n * B, replace = TRUE), nrow = n)
  } else if (bootstrap_method == "stationary") {
    if (is.null(block_length)) stop("No average block length provided!")
    idx <- stationary_bootstrap_indices(n = n, avg_block_size = block_length, B = B)
  } else {
    stop("Not a valid bootstrap method")
  }

  # Estimate the model on the bootstraped data
  if (methods::is(fit, "esreg")) {
    # Use the one_shot estimation approach for speed
    b <- apply(idx, 2, function(id) {
      fitb <- esreg(fit$y[id] ~ fit$x[id, -1],
                    alpha = fit$alpha, g1 = fit$g1, g2 = fit$g2, method = 'one_shot')
      fitb$par
    })
  } else if (methods::is(fit, "esreg_twostep")) {
    b <- apply(idx, 2, function(id) {
      fitb <- esreg_twostep(fit$y[id] ~ fit$x[id, -1], alpha = fit$alpha)
      fitb$par
    })
  }

  # Compute the covariance
  cov <- stats::cov(t(b))

  cov
}
