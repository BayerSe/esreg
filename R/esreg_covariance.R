#' Bandwidth
#'
#' Bandwidth for the estimation of the density quantile function
#'
#' @param n sample size
#' @param alpha quantile index
#' @param tau alpha level of the test
#' @param type Bofinger, Chamberlain or Hall-Sheather
#' @references Bofinger (1975), Chamberlain (1994), Hall and Sheather(1988)
#' @export
bandwidth <- function(n, alpha, tau = 0.05, type = "Hall-Sheather") {
  if (type == "Bofinger") {
    n^(-1/5) * ((9/2 * stats::dnorm(stats::qnorm(alpha))^4)/(2 * stats::qnorm(alpha)^2 + 1)^2)^(1/5)
  } else if (type == "Chamberlain") {
    stats::qnorm(1 - alpha/2) * sqrt(tau * (1 - tau)/n)
  } else if (type == "Hall-Sheather") {
    n^(-1/3) * stats::qnorm(1 - tau/2)^(2/3) * ((3/2 * stats::dnorm(stats::qnorm(alpha))^2)/(2 * stats::qnorm(alpha)^2 + 1))^(1/3)
  }
}


#' Conditional truncated variance
#'
#' Estimate the conditional truncated variance
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param approach scl_N or scl_t
#'
#' @export
conditional_truncated_variance <- function(y, x, approach) {

  # Some variables
  z <- x  # Variables relevant for the standard deviation
  n <- nrow(x)
  k <- ncol(x)
  l <- ncol(z)
  scaling <- strsplit(approach, "_")[[1]][2]

  # Mean / standard deviation specifications
  mu_fun <- function(x, b) as.numeric(x %*% b)
  sigma_fun <- function(x, b) as.numeric(x %*% b)

  # Normal log-likelihood
  f <- function(b, y, x, z) {
    mu <- mu_fun(x, b = b[1:k])
    sigma <- sigma_fun(z, b = b[(k + 1):(k + l)])
    if (all(sigma > 0)) {
      stats::dnorm(y, mean = mu, sd = sigma, log = TRUE)
    } else {
      rep(NA, n)
    }
  }

  # Starting values (and ensure positive fitted standard deviations)
  fit1 <- stats::lm(y ~ x - 1)
  fit2 <- stats::lm(abs(fit1$residuals) ~ z - 1)
  fit2$coefficients[1] <- fit2$coefficients[1] - min(0.001, min(fit2$fitted.values))
  b0 <- c(fit1$coefficients, fit2$coefficients)

  # Optimize the model
  fit <- NULL
  for (method in c("BFGS", "BHHH", "NM")) {
    fit <- try(maxLik::maxLik(f, start = b0, y = y, x = x, z = z, method = method,
                              control = list(iterlim = 1000)), silent = TRUE)
    if (inherits(fit, "maxLik"))
      if (fit$code == 0)
        break
  }

  # Fitted values
  b <- fit$estimate
  mu <- mu_fun(x, b = b[1:k])
  sigma <- sigma_fun(z, b = b[(k + 1):(k + l)])

  # Refinements and Truncation
  if (scaling == "N") {
    beta <- -mu / sigma
    beta[beta < -30] <- -30
    cv <- sigma^2 * (1 - beta * stats::dnorm(beta)/stats::pnorm(beta) - (stats::dnorm(beta)/stats::pnorm(beta))^2)
  } else if (scaling == "t") {
    # Student-t log likelihood
    f <- function(b, y, x, z) {
      mu <- mu_fun(x, b = b[1:k])
      sigma <- sigma_fun(z, b = b[(k + 1):(k + l)])
      df <- utils::tail(b, 1)
      if (all(sigma > 0) & (df > 2.1)) {
        lgamma((df + 1)/2) - lgamma(df/2) - log(sqrt(pi * df) * sigma) - (df + 1)/2 * log((1 + 1/df * ((y - mu)/sigma)^2))
      } else {
        rep(NA, n)
      }
    }

    # Inital estimate of the degrees of freedom
    r <- (y - mu) / sigma
    df <- stats::optimise(function(df, x = r, mu = 0, sigma = 1) {
      sigma <- sigma * sqrt((df - 2)/df)  # Correct the standard deviation
      -sum(lgamma((df + 1)/2) - lgamma(df/2) - log(sqrt(pi * df) * sigma) - (df + 1)/2 * log((1 + 1/df * ((x - mu)/sigma)^2)))
    }, interval = c(2.1, 1e+06))$minimum

    # Estimate the full model
    fit <- NULL
    for (method in c("BFGS", "BHHH", "NM")) {
      fit <- try(maxLik::maxLik(f, start = c(b, df), y = y, x = x, z = z, method = method,
                                control = list(iterlim = 1000)), silent = TRUE)
      if (inherits(fit, "maxLik"))
        if (fit$code == 0)
          break
    }

    # Fitted values
    b <- fit$estimate
    mu <- mu_fun(x, b = b[1:k])
    sigma <- sigma_fun(z, b = b[(k + 1):(k + l)])
    df <- utils::tail(b, 1)
    sigma <- sigma * sqrt(df/(df - 2))  # Transform back to the standard deviation

    # See Ho et al. (2012): Some results on the truncated multivariate t distribution
    t1 <- (df - 2) / df
    beta <- -mu / (sqrt(t1) * sigma)
    beta[beta < -30] <- -30
    if (df < 300) {
      k <- gamma((df+1)/2) / gamma(df/2) / sqrt(df*pi) / stats::pt(beta, df=df)
    } else {
      k <- sqrt(df/2) / sqrt(df*pi) / stats::pt(beta, df=df)
    }
    m1 <- k * df / (df - 1) * (-(1+beta^2/df)^(-(df-1)/2))
    m2 <- (df - 1) / t1 * (stats::pt(beta*sqrt(t1), df=df-2) / stats::pt(beta, df=df)) - df
    cv <- t1 * sigma^2 * (m2 - m1^2)
  }

  # Return
  cv
}


#' Asymptotic covariance
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

  # Check the inputs
  if (!(sparsity %in% c("iid", "nid")))
    stop("sparsity can be iid or nid")
  if (!(cond_var %in% c('ind', 'scl_N', 'scl_t')))
    stop('cond_var can be ind, scl_N or scl_t')
  if (!(bandwidth_type %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
    stop("bandwidth_type can be Bofinger, Chamberlain or Hall-Sheather")

  # Extract some elements from the fit object and compute some things
  y <- fit$y
  x <- fit$x
  alpha <- fit$alpha
  par_q <- fit$par_q
  par_e <- fit$par_e
  xq <- as.numeric(x %*% par_q)
  xe <- as.numeric(x %*% par_e)
  u <- as.numeric(y - xq)
  eps <- .Machine$double.eps^(2/3)
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

  # Transform the data and coefficients
  if (fit$g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
    par_q[1] <- par_q[1] - max_y
    par_e[1] <- par_e[1] - max_y
  }

  # Density quantile function
  if (sparsity == "iid") {
    # Koenker (1994)
    h <- max(k + 1, ceiling(n * bandwidth(n = n, alpha = alpha, type = bandwidth_type)))
    ir <- (k + 1):(h + k + 1)
    ord.resid <- sort(u[order(abs(u))][ir])
    xt <- ir/(n - k)
    dens <- 1/suppressWarnings(quantreg::rq(ord.resid ~ xt)$coef[2])
    dens <- rep(as.numeric(dens), n)
  } else if (sparsity == "nid") {
    # Hendricks and Koenker (1992)
    h <- bandwidth(n = n, alpha = alpha, type = bandwidth_type)
    bu <- quantreg::rq(y ~ x - 1, tau = alpha + h)$coef
    bl <- quantreg::rq(y ~ x - 1, tau = alpha - h)$coef
    dens <- pmax(0, 2 * h/(x %*% (bu - bl) - eps))
  }

  # Truncated conditional variance
  if (cond_var == "ind") {
    cv <- rep(stats::var(u[u <= 0]), n)
  } else {
    cv <- tryCatch({
      cv <- conditional_truncated_variance(y = u, x = x, approach = cond_var)
      if (any(is.na(cv) | any(!is.finite(cv))))
        stop() else cv
    }, error = function(e) {
      warning(paste0("Can not fit ", cond_var, ", now using ind estimator!"))
      rep(stats::var(u[u <= 0]), n)
    })
  }

  # Evaluate G1 / G2 functions
  G1_prime_xq <- G_vec(z = xq, g = 'G1_prime', type = fit$g1)
  G2_xe <- G_vec(z = xe, g = 'G2', type = fit$g2)
  G2_prime_xe <- G_vec(z = xe, g = 'G2_prime', type = fit$g2)

  # Compute the covariance
  loop <- loop_esreg_covariance(x = x, xq = xq, xe = xe, alpha = alpha,
                                G1_prime_xq = G1_prime_xq,
                                G2_xe = G2_xe, G2_prime_xe = G2_prime_xe,
                                density = dens, conditional_variance = cv)

  # Check if we can compute the covariance matrix, if not return NA matrix
  cov_matrix <- suppressWarnings(try(solve(loop$lambda) %*% loop$C %*% solve(loop$lambda)))
  if (class(cov_matrix) == "try-error") {
    matrix(NA, 2 * k, 2 * k)
  } else {
    rownames(cov_matrix) <- colnames(cov_matrix) <- names(fit$par)
    cov_matrix
  }
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
  # Use the one_shot estimation approach for speed
  b <- apply(idx, 2, function(id) {
    fitb <- esreg(fit$y[id] ~ fit$x[id, -1],
                  alpha = fit$alpha, g1 = fit$g1, g2 = fit$g2, method = 'one_shot')
    fitb$par
  })

  # Compute the covariance
  cov <- stats::cov(b)

  cov
}
