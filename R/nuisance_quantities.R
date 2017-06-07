#' @title Denstiy Quantile Function
#' @description Estimate the density quantile function
#' @param y Vector of dependent data
#' @param x Matrix of covariates
#' @param u Quantile residuals
#' @param alpha Quantile of interest
#' @param sparsity iid or ind
#' @param bandwidth_type Bofinger, Chamberlain or Hall-Sheather
#' @references Bofinger (1975), Chamberlain (1994), Hall and Sheather(1988)
#' @export
density_quantile_function <- function(y, x, u, alpha, sparsity, bandwidth_type) {
  if (!("quantreg" %in% rownames(utils::installed.packages()))) {
    stop("quantreg needed for this function to work. Please install it.")
  }

  n <- nrow(x)
  k <- ncol(x)
  eps <- .Machine$double.eps^(2/3)

  # Obtain the bandwidth
  if (bandwidth_type == "Bofinger") {
    bandwidth <- n^(-1/5) * ((9/2 * stats::dnorm(stats::qnorm(alpha))^4)/(2 * stats::qnorm(alpha)^2 + 1)^2)^(1/5)
  } else if (bandwidth_type == "Chamberlain") {
    tau <- 0.05
    bandwidth <- stats::qnorm(1 - alpha/2) * sqrt(tau * (1 - tau)/n)
  } else if (bandwidth_type == "Hall-Sheather") {
    tau <- 0.05
    bandwidth <- n^(-1/3) * stats::qnorm(1 - tau/2)^(2/3) * ((3/2 * stats::dnorm(stats::qnorm(alpha))^2)/(2 * stats::qnorm(alpha)^2 + 1))^(1/3)
  } else {
    stop("Not a valid bandwith method!")
  }

  # Compute the density
  if (sparsity == "iid") {
    # Koenker (1994)
    h <- max(k + 1, ceiling(n * bandwidth))
    ir <- (k + 1):(h + k + 1)
    ord.resid <- sort(u[order(abs(u))][ir])
    xt <- ir/(n - k)
    density <- 1/suppressWarnings(quantreg::rq(ord.resid ~ xt)$coef[2])
    density <- rep(as.numeric(density), n)
  } else if (sparsity == "nid") {
    # Hendricks and Koenker (1992)
    bu <- quantreg::rq(y ~ x - 1, tau = alpha + bandwidth)$coef
    bl <- quantreg::rq(y ~ x - 1, tau = alpha - bandwidth)$coef
    density <- pmax(0, 2 * bandwidth/(x %*% (bu - bl) - eps))
  } else {
    stop("Not a valid density quantile function estimator!")
  }

  density
}


#' @title Conditional truncated variance
#' @description Estimate the conditional truncated variance of y given x
#' and given y <= 0.
#' @param y Vector of dependent data
#' @param x Matrix of covariates including the intercept
#' @param approach ind, scl_N, scl_t or scl_st
#' @param approximate Concerns the scl_st estimator: approximate the time consuming
#' integration via splines.
#' @export
conditional_truncated_variance <- function(y, x, approach, approximate=TRUE) {

  if (sum(y <= 0) <= 2) {
    stop("Not enough negative quantile residuals!")
  }

  if (approach == "ind") {
    cv <- rep(stats::var(y[y <= 0]), length(y))
  } else {
    cv <- tryCatch({
      if (approach == "scl_N") {
        cv <- .conditional_truncated_variance_N(y=y, x=x)
      } else if (approach == "scl_t") {
        cv <- .conditional_truncated_variance_t(y=y, x=x)
      } else if (approach == "scl_st") {
        cv <- .conditional_truncated_variance_st(y=y, x=x, approximate=TRUE)
      } else {
        stop("Not a valid estimator!")
      }
      if (any(is.na(cv) | any(!is.finite(cv)))) stop() else cv
    }, error = function(e) {
      warning(paste0("Can not fit the ", approach, " estimator, switching to the ind approach!"))
      rep(stats::var(y[y <= 0]), length(y))
    })
  }

  cv
}

.conditional_truncated_variance_N <- function(y, x, return_estimates = FALSE) {
  # Obtain starting values and ensure positive fitted standard deviations
  fit1 <- stats::lm(y ~ x - 1)
  fit2 <- stats::lm(abs(fit1$residuals) ~ x - 1)
  fit2$coefficients[1] <- fit2$coefficients[1] - min(0.001, min(fit2$fitted.values))
  b0 <- c(fit1$coefficients, fit2$coefficients)

  # Estimate the full model
  ll <- function(par, y, x) {
    k <- ncol(x)
    mu <- as.numeric(x %*% par[1:k])
    sigma <- as.numeric(x %*% par[(k+1):(2*k)])
    if (all(sigma > 0)) {
      sum(-stats::dnorm(x=y, mean=mu, sd=sigma, log=TRUE))
    } else {
      NA
    }
  }
  b <- .mle_model(b0=b0, y=y, x=x, ll=ll)

  # Estimated means and standard deviations
  k <- ncol(x)
  mu <- as.numeric(x %*% b[1:k])
  sigma <- as.numeric(x %*% b[(k+1):(2*k)])

  # Return results
  if (return_estimates) {
    # Return the estimates as input for the other estimators
    list(b=b, mu=mu, sigma=sigma)
  } else {
    # Truncated variance
    beta <- -mu / sigma
    beta[beta < -30] <- -30
    cv <- sigma^2 * (1 - beta * stats::dnorm(beta)/stats::pnorm(beta) - (stats::dnorm(beta)/stats::pnorm(beta))^2)

    cv
  }
}

.conditional_truncated_variance_t <- function(y, x) {
  # Get the estimates under the normality assumption as starting values
  fit0 <- .conditional_truncated_variance_N(y, x, return_estimates = TRUE)
  b0 <- fit0$b
  mu <- fit0$mu
  sigma <- fit0$sigma

  # Inital estimate of the degrees of freedom
  df <- stats::optimise(function(df, x = scale((y - mu) / sigma)) {
    sigma <- sqrt((df-2)/df)
    -sum(lgamma((df+1)/2) - lgamma(df/2) - log(sqrt(pi*df) * sigma) - (df+1)/2 * log(1 + 1/df * (x/sigma)^2))
  }, interval = c(2.1, 1e+06))$minimum


  # Estimate the full model
  ll <- function(par, y, x) {
    k <- ncol(x)
    mu <- as.numeric(x %*% par[1:k])
    sigma <- as.numeric(x %*% par[(k+1):(2*k)])
    df <- par[2*k+1]
    if (all(sigma > 0) & (df > 2.1)) {
      sigma <- sigma * sqrt((df-2)/df)
      -sum(lgamma((df+1)/2) - lgamma(df/2) - log(sqrt(pi*df) * sigma) - (df+1)/2 * log((1 + 1/df * ((y - mu)/sigma)^2)))
    } else {
      NA
    }
  }
  b <- .mle_model(b0=c(b0, df), y=y, x=x, ll=ll)

  # Estimated means, standard deviations and the degrees of freedom
  k <- ncol(x)
  mu <- as.numeric(x %*% b[1:k])
  sigma <- as.numeric(x %*% b[(k+1):(2*k)])
  df <- b[2*k+1]

  # Truncated variance, see Ho et al. (2012): Some results on the truncated multivariate t distribution
  t1 <- (df-2) / df
  beta <- -mu / (sqrt(t1) * sigma)
  beta[beta < -30] <- -30
  if (df < 300) {
    k <- gamma((df+1)/2) / gamma(df/2) / sqrt(df*pi) / stats::pt(beta, df=df)
  } else {
    k <- sqrt(df/2) / sqrt(df*pi) / stats::pt(beta, df=df)
  }
  m1 <- k * df / (df-1) * (-(1+beta^2/df)^(-(df-1)/2))
  m2 <- (df - 1) / t1 * (stats::pt(beta*sqrt(t1), df=df-2) / stats::pt(beta, df=df)) - df
  cv <- t1 * sigma^2 * (m2 - m1^2)

  cv
}

.conditional_truncated_variance_st <- function(y, x, approximate=TRUE) {
  # Get the estimates under the normality assumption as starting values
  fit0 <- .conditional_truncated_variance_N(y, x, return_estimates = TRUE)
  b0 <- fit0$b
  mu <- fit0$mu
  sigma <- fit0$sigma

  # Inital estimate of the degrees of freedom and the skewness
  fit_init <- stats::optim(par=c(5, 1), method="Nelder-Mead",
                           fn=function(b) -sum(.sdst(x=scale((y - mu) / sigma),
                                                     df=b[1], skew=b[2], log=TRUE)))
  df <- fit_init$par[1]
  skew <- fit_init$par[2]

  # Estimate the full model
  ll <- function(par, y, x) {
    k <- ncol(x)
    mu <- as.numeric(x %*% par[1:k])
    sigma <- as.numeric(x %*% par[(k+1):(2*k)])
    df <- par[2*k+1]
    skew <- par[2*k+2]
    if (all(sigma > 0) & (df > 2.1) & (skew > 0.01)) {
      -sum(.sdst(x=(y-mu)/sigma, df=df, skew=skew, log=TRUE) - log(sigma))
    } else {
      NA
    }
  }
  b <- .mle_model(b0=c(b0, df, skew), y=y, x=x, ll=ll)

  # Estimated means, standard deviations, degrees of freedom and the skewness
  k <- ncol(x)
  mu <- as.numeric(x %*% b[1:k])
  sigma <- as.numeric(x %*% b[(k+1):(2*k)])
  df <- b[2*k+1]
  skew <- b[2*k+2]

  # Truncated variance
  td <- function(x, b, df, skew) .sdst(x=x, df=df, skew=skew, log=FALSE) / .spst(x=b, df=df, skew=skew)
  beta <- -mu / sigma
  beta[beta < -30] <- -30
  if (approximate) {
    b_approx <- seq(min(beta), max(beta), length.out = 50)
    cv_approx <- sapply(b_approx, function(b) {
      m1 <- stats::integrate(function(x) x   * td(x, b, df, skew), lower = -Inf, upper = b)$value
      m2 <- stats::integrate(function(x) x^2 * td(x, b, df, skew), lower = -Inf, upper = b)$value
      m2 - m1^2
    })
    cv <- sigma^2 * stats::spline(x=b_approx, y=cv_approx, xout=beta)$y
  } else {
    cv <- sigma^2 * sapply(beta, function(b) {
      m1 <- stats::integrate(function(x) x   * td(x, b, df, skew), lower = -Inf, upper = b)$value
      m2 <- stats::integrate(function(x) x^2 * td(x, b, df, skew), lower = -Inf, upper = b)$value
      m2 - m1^2
    })
  }

  cv
}

# Estimate the models via maximum likelihood estimation
.mle_model <- function(b0, y, x, ll) {
  fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="BFGS"), silent=TRUE)
  if(inherits(fit, "try-error") || (fit$convergence != 0)) {
    fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="Nelder-Mead",
                            control=list(maxit=10000)), silent=TRUE)
  }
  if(inherits(fit, "try-error") || (fit$convergence != 0)) {
    stop("Cannot fit the model!")
  }
  fit$par
}

# Standardized Student-t distribution (mean 0, variance 1)
.dst <- function(x, df, log=FALSE) {
  sigma <- sqrt((df - 2)/df)
  l <- lgamma((df+1)/2) - lgamma(df/2) - log(sqrt(pi*df) * sigma) -
    (df+1)/2 * log(1 + 1/df * (x/sigma)^2)
  if (log) { l } else { exp(l) }
}
.pst <- function(x, df) {
  stats::pt(x / sqrt((df-2)/df), df)
}

# Standardized skewed Student-t distribution using the Fernarndez-Steel approach
.sdst <- function(x, df, skew, log=FALSE) {
  if ((df < 2.1) | (skew < 0.1)) return(NA)
  m1 <- 2 * sqrt(df-2) / (df-1) / beta(0.5, 0.5 * df)
  m <- m1 * (skew - 1/skew)
  s <- sqrt((1 - m1^2) * (skew^2 + 1/skew^2) + 2*m1^2 - 1)
  z <- m + s*x
  l <- log(2*s / (skew + 1/skew)) +
    (.dst(z/skew, df=df, log=TRUE) * (x >= -m/s) + .dst(z*skew, df=df, log=TRUE) * (x < -m/s))
  if (log) { l } else { exp(l) }
}
.spst <- function(x, df, skew) {
  m1 <- 2 * sqrt(df-2) / (df-1) / beta(0.5, 0.5 * df)
  m <- m1 * (skew - 1/skew)
  s <- sqrt((1 - m1^2) * (skew^2 + 1/skew^2) + 2*m1^2 - 1)
  z <- m + s*x

  p1 <- (2 / (skew + 1/skew) * (skew*.pst(z / skew, df) + 1/skew) - 1) * (x >= -m/s)
  p2 <- (2 / skew / (skew + 1/skew) * .pst(z * skew, df)) * (x < -m/s)

  p1 + p2
}
