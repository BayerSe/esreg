#' Joint (VaR, ES) Regression
#'
#' @param formula y ~ x1 + x2 + ...
#' @param data data.frame that stores y and x. Extracted from the enviroment if missing.
#' @param alpha Quantile index
#' @param g1 1, 2 (see \link{G1_fun}, \link{G1_prime_fun})
#' @param g2 1, 2, 3, 4, 5 (see \link{G2_curly_fun}, \link{G2_fun}, \link{G2_prime_fun})
#' @param b0 Starting values; if NULL they are obtained from two additional quantile regressions
#' @param target The functions to be optimized: either the loss (rho) or the identification function (psi).
#' @param method one_shot, random_restart or gensa
#' @param random_restart_ctrl list: N number of random starting points; M optimize over the M best; sd std dev of the random component
#' @param gensa_ctrl Parameters to be passed to GenSA
#' @export
esreg <- function(formula, data, alpha, g1 = 2L, g2 = 1L, b0 = NULL, target = "rho", method = "random_restart",
                  random_restart_ctrl = list(M = 10, N = 1000, sd = sqrt(0.1)),
                  gensa_ctrl = list(max.call = 1e+06, max.time = 10, box = 10)) {

  # Check the inputs
  if (!(g1 %in% c(1, 2)))
    stop("G1 can be 1 or 2.")
  if (!(g2 %in% c(1, 2, 3, 4, 5)))
    stop("G2 can be 1, 2, 3, 4 or 5.")
  if ((alpha < 0) | (1 < alpha))
    stop("alpha not in (0, 1)")
  if (!(target %in% c("rho", "psi")))
    stop("target can be rho or psi")
  if (!(method %in% c("one_shot", "random_restart", "gensa")))
    stop("method can be one_shot, random_restart or gensa")

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

  # Set some variables
  k <- ncol(x)
  optim_ctrol <- list(maxit = 10000, reltol = .Machine$double.eps^(1/2))

  # Transform the data
  if (g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
  }

  # Find starting values
  if (is.null(b0)) {
    if (k == 1) {
      q <- stats::quantile(y, p = alpha, type = 8)
      b0 <- as.numeric(c(q, mean(y[y <= q])))
    } else {
      e <- -stats::dnorm(stats::qnorm(alpha))/alpha  # Match quantile and expected shortfall
      alpha_tilde <- stats::uniroot(function(x) stats::qnorm(x) - e, c(0, alpha))$root
      b0 <- as.vector(quantreg::rq(y ~ x - 1, tau = c(alpha, alpha_tilde))$coef)
    }
  } else {
    if (length(b0) != 2 * k)
      stop("Length of b0 does not match to k!")
  }

  # Split the starting vector
  b0_q <- b0[1:k]
  b0_e <- b0[(k + 1):(2 * k)]

  # Target functions
  if (target == "rho") {
    fun1 <- function(b) suppressWarnings(esr_rho_lp(b = c(b0_q, b), y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
    fun2 <- function(b) suppressWarnings(esr_rho_lp(b = b, y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
  } else if (target == "psi") {
    fun1 <- function(b) suppressWarnings(esr_psi_lp(b = c(b0_q, b), y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
    fun2 <- function(b) suppressWarnings(esr_psi_lp(b = b, y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
  }

  # Optimize the model
  if (k == 1) {
    fit <- stats::optim(par = b0, fn = fun2, method = "Nelder-Mead", control = optim_ctrol)
    method <- "direct_optimization"
  } else {
    if (method == "one_shot") {
      # Refine the shortfall starting values
      fit0 <- stats::optim(par = b0_e, fn = fun1, method = "Nelder-Mead",
                           control = optim_ctrol)
      b0 <- c(b0_q, fit0$par)
      fit <- stats::optim(par = b0, fn = fun2, method = "Nelder-Mead",
                          control = optim_ctrol)
    } else if (method == "random_restart") {
      # Evaluate the N random starting points
      noise <- matrix(stats::rnorm(2*k * random_restart_ctrl$N, sd = random_restart_ctrl$sd),
                      nrow = random_restart_ctrl$N)
      rand_start <- sweep(noise, MARGIN = 2, b0, "+")
      rand_start <- cbind(apply(rand_start, 1, fun2), rand_start)

      # Optimize over the M best and the original starting value
      start_values <- rbind(rand_start[order(rand_start[, 1])[1:random_restart_ctrl$M], -1], b0)
      fits <- apply(start_values, 1, function(b)
        stats::optim(par = b, fn = fun2, method = "Nelder-Mead", control = optim_ctrol))

      # Find the best fit
      fit <- fits[[which.min(sapply(fits, "[[", "value"))]]
    } else if (method == "gensa") {
      if (!("GenSA" %in% rownames(utils::installed.packages()))) {
        stop("GenSA needed for this function to work. Please install it.")
      }
      set.seed(1)
      fit <- GenSA::GenSA(par = b0, fn = fun2,
                          lower = b0 - rep(gensa_ctrl$box, 2 * k),
                          upper = b0 + rep(gensa_ctrl$box, 2 * k),
                          control = gensa_ctrl[!(names(gensa_ctrl) %in% c("box"))])
    }
  }

  # Set names of the parameters
  b <- fit$par
  names(b) <- c(paste0("bq_", 1:k - 1), paste0("be_", 1:k - 1))

  # Undo the transformation
  if (g2 %in% c(1, 2, 3)) {
    y <- y + max_y
    b[1] <- b[1] + max_y
    b[k + 1] <- b[k + 1] + max_y
    b0[1] <- b0[1] + max_y
    b0[k + 1] <- b0[k + 1] + max_y
  }

  # Return results
  structure(list(call = cl, target = target, method = method, g1 = g1, g2 = g2,
                 alpha = alpha, y = y, x = x, b0 = b0,
                 par = b, par_q = b[1:k], par_e = b[(k + 1):(2 * k)],
                 value = fit$value, time = Sys.time() - t0),
            class = "esreg")
}

#' @export
print.esreg <- function(x, ...) {
  cat("Method:", x$method, "\n")
  cat("Target:", x$target, "\n")
  cat(.G_function_names(x$g1, x$g2)[1], "\n")
  cat(.G_function_names(x$g1, x$g2)[2], "\n")
  cat("Value: ", sprintf("%.9f", x$value), "\n")
  cat("alpha: ", sprintf("%.3f", x$alpha), "\n")
  cat("Time:  ", sprintf("%.3f", x$time), "\n\n")
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  cat("Estimates:\n")
  cat(sprintf("% 0.4f", x$par_q), "\n")
  cat(sprintf("% 0.4f", x$par_e))
}

#' @export
coef.esreg <- function(object, ...) {
  object$par
}

#' @export
fitted.esreg <- function(object, ...) {
  cbind(object$x %*% object$par_q, object$x %*% object$par_e)
}

#' Estimated covariance of the joint (VaR, ES) estimator
#'
#' Estimate the variance-covariance matrix of the joint (VaR, ES) estimator
#' either using the asymptotic formulas or using the bootstrap.
#'
#' @param sparsity Sparsity estimator
#' \itemize{
#'   \item iid - Piecewise linear interpolation of the distribution
#'   \item nid - Hendricks and Koenker sandwich
#' }
#' @param cond_var Conditional truncated variance estimator
#' \itemize{
#'   \item ind Variance over all negative residuals
#'   \item scl_N Scaling with the Normal distribution
#'   \item scl_t Scaling with the t-distribution
#' }
#' @param bandwidth_type Bofinger, Chamberlain or Hall-Sheather
#' @param bootstrap_method
#' \itemize{
#'   \item NULL asymptotic estimator
#'   \item iid
#'   \item stationary Politis & Romano (1994)
#' }
#' @param B Number of bootstrap iterations
#' @param block_length Average block length for the stationary bootstrap
#' @export
vcov.esreg <- function(object, sparsity = "iid", cond_var = "ind", bandwidth_type = "Hall-Sheather",
                       bootstrap_method = NULL, B = 1000, block_length = NULL, ...) {
  fit <- object

  if(is.null(bootstrap_method)) {
    if(!(sparsity %in% c("iid", "nid")))
      stop("sparsity can be iid or nid")
    if(!(cond_var %in% c('ind', 'scl_N', 'scl_t')))
      stop('cond_var can be ind, scl_N or scl_t')
    if(!(bandwidth_type %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
      stop("bandwidth_type can be Bofinger, Chamberlain or Hall-Sheather")

    # Extract some elements
    y <- fit$y
    x <- fit$x
    n <- nrow(x)
    k <- ncol(x)
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

    # Check the methods in case of sample quantile / es
    if ((k == 1) & sparsity != "iid") {
      warning("Changed sparsity estimation to iid!")
      sparsity <- "iid"
    }
    if ((k == 1) & cond_var != "ind") {
      warning("Changed conditional truncated variance estimation to nid!")
      cond_var <- "ind"
    }

    # Density quantile function
    dens <- density_quantile_function(y = y, x = x, u = u, alpha = fit$alpha,
                                      sparsity = sparsity, bandwidth_type = bandwidth_type)

    # Truncated conditional variance
    cv <- conditional_truncated_variance(y = y, x = x, u = u, approach = cond_var)

    # Evaluate G1 / G2 functions
    G1_prime_xq <- G_vec(z = xq, g = "G1_prime", type = fit$g1)
    G2_xe <- G_vec(z = xe, g = "G2", type = fit$g2)
    G2_prime_xe <- G_vec(z = xe, g = "G2_prime", type = fit$g2)

    # Compute the covariance matrix
    cov <- l_esreg_covariance(x = x, xq = xq, xe = xe, alpha = fit$alpha,
                              G1_prime_xq = G1_prime_xq,
                              G2_xe = G2_xe, G2_prime_xe = G2_prime_xe,
                              density = dens, conditional_variance = cv)
  } else {
    if (!(bootstrap_method %in% c(NULL, "iid", "stationary")))
      stop("bootstrap_method can be NULL, iid or stationary")
    if (B < 1000)
      warning("The number of bootstrap iterations is small!")

    # Draw the bootstrap indices
    n <- length(fit$y)
    set.seed(1)
    if (bootstrap_method == "iid") {
      idx <- matrix(sample(1:n, size = n * B, replace = TRUE), nrow = n)
    } else if (bootstrap_method == "stationary") {
      if (is.null(block_length))
        stop("No average block length provided!")
      idx <- stationary_bootstrap_indices(n = n, avg_block_size = block_length, B = B)
    }

    # Use the one_shot estimation approach for speed
    b <- apply(idx, 2, function(id) {
      fitb <- esreg(fit$y[id] ~ fit$x[id, -1],
                    alpha = fit$alpha, g1 = fit$g1, g2 = fit$g2, method = 'one_shot')
      fitb$par
    })

    # Compute the covariance
    cov <- stats::cov(t(b))
  }

  # Set the names
  rownames(cov) <- colnames(cov) <- names(coef(fit))

  # Return the estimated covariance
  cov
}
