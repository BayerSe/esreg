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
#' @description Estimate the conditional truncated variance of y given x when y is truncated at 0.
#' @param y Vector of dependent data
#' @param x Matrix of covariates
#' @param approach ind, scl_N or scl_t
#' @export
conditional_truncated_variance <- function(y, x, approach) {
  if (sum(y <= 0) <= 2) {
    stop("Not enough negative quantile residuals!")
  }

  if (approach == "ind") {
    cv <- rep(stats::var(y[y <= 0]), length(y))
  } else if (approach %in% c("scl_N", "scl_t")) {
    if (!("maxLik" %in% rownames(utils::installed.packages()))) {
      stop("maxLik needed for this function to work. Please install it.")
    }
    cv <- tryCatch({
      # Store some variables
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

      if (any(is.na(cv) | any(!is.finite(cv)))) stop() else cv
    }, error = function(e) {
      warning(paste0("Can not fit the ", approach, " estimator, I'm switching to the ind approach!"))
      rep(stats::var(y[y <= 0]), length(y))
    })
  } else {
    stop("Not a valid estimator!")
  }

  cv
}
