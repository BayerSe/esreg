#' Two Step (VaR, ES) Regression
#'
#' Estimates the expected shortfall in two steps.
#'
#' This estimator is much faster than the
#' one-step estimator \link{esreg}. Its estimates are, however,
#' less precise and therefore \link{esreg} is geneally the preferred estimator.
#'
#' @param formula y ~ x1 + x2 + ...
#' @param data data.frame that stores y and x. Extracted from the enviroment if missing.
#' @param alpha Quantile index
#' @export
esreg_twostep <- function(formula, data, alpha) {

  # Start the timer
  t0 <- Sys.time()

  # Extract the formula
  if (missing(data))
    data <- environment(formula)
  cl <- match.call()
  mf <- stats::model.frame(formula = formula, data = data)
  x <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y <- stats::model.response(mf)

  # Check the data
  if (any(is.na(y)) | any(is.na(x)))
    stop("Data contains NAs!")
  if (!(all(is.finite(y)) & all(is.finite(x))))
    stop("Not all values are finite!")

  # First stage: quantile regression
  fit_rq <- quantreg::rq(y ~ x-1, tau = alpha)
  q <- fit_rq$fitted.values
  b_q <- fit_rq$coefficients

  # Second stage: weighted least squares regression
  fit <- stats::lm(y ~ x-1, weights = (y <= q) * 1)
  b_e <- fit$coefficients

  # Name the coefficents
  names(b_q) <- paste0("bq_", 1:length(b_q) - 1)
  names(b_e) <- paste0("be_", 1:length(b_e) - 1)

  # Return results
  structure(list(call = cl, alpha = alpha, y = y, x = x,
                 par = c(b_q, b_e), par_q = b_q, par_e = b_e,
                 time = Sys.time() - t0),
            class = "esreg_twostep")
}

#' @export
print.esreg_twostep <- function(x, ...) {
  cat("alpha: ", sprintf("%.3f", x$alpha), "\n")
  cat("Time:  ", sprintf("%.3f", x$time), "\n\n")
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  cat("Estimates:\n")
  cat(sprintf("% 0.4f", x$par_q), "\n")
  cat(sprintf("% 0.4f", x$par_e))
}

#' @export
coef.esreg_twostep <- function(object, ...) {
  object$par
}

#' @export
fitted.esreg_twostep <- function(object, ...) {
  cbind(object$x %*% object$par_q, object$x %*% object$par_e)
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
