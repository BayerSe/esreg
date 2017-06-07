#' @title Joint (VaR, ES) Regression
#' @description  Estimates a joint linear regression model for the pair (VaR, ES):
#' \deqn{Q_\alpha(Y | X) = X'\beta_q}
#' \deqn{ES_\alpha(Y | X) = X'\beta_e}
#' @param formula Forumula object, e.g.: y ~ x1 + x2 + ...
#' @param data data.frame that holds the variables. Can be missing.
#' @param alpha Quantile of interest
#' @param g1 1, [2] (see \link{G1_fun}, \link{G1_prime_fun})
#' @param g2 [1], 2, 3, 4, 5 (see \link{G2_curly_fun}, \link{G2_fun}, \link{G2_prime_fun})
#' @param b0 Starting values for the optimization; if NULL they are obtained from two additional quantile regressions
#' @param target The functions to be optimized: either the loss [rho] or the identification function (psi). The rho function is strongly recommended.
#' @param method [random_restart] or gensa
#' @param shift_data If g2 is 1, 2 or 3, we can either estimate the model without or with
#' shifting of the Y variable. We either risk biased estimates (no shifting) or slightly different estimates due
#' to the changed loss function (shifting). Defaults to shifting to avoid biased estimates.
#' @param random_restart_ctrl list: N [1000] number of random starting points; M [10] optimize over the M best; sd [sqrt(0.1)] std dev of the random component
#' @param gensa_ctrl Parameters to be passed to the GenSA opzimizer
#' @return An esreg object
#' @seealso \code{\link{vcov.esreg}} for the covariance estimation and
#' \code{\link{summary.esreg}} for a summary of the regression results
#' @examples
#' # Simulate data
#' set.seed(1)
#' x <- rchisq(1000, df = 1)
#' y <- -x + (1 + 0.1 * x) * rnorm(1000)
#'
#' # Estimate the model and the covariance
#' fit <- esreg(y ~ x, alpha = 0.025)
#' summary(object = fit, sparsity = "nid", cond_var = "scl_t")
#' @references \href{https://arxiv.org/abs/1704.02213}{A Joint Quantile and Expected Shortfall Regression Framework}
#' @export
esreg <- function(formula, data, alpha, g1 = 2L, g2 = 1L, b0 = NULL, target = "rho", method = "random_restart",
                  shift_data = TRUE, random_restart_ctrl = list(M = 10, N = 1000, sd = sqrt(0.1)),
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
  if (!(method %in% c("random_restart", "gensa")))
    stop("method can be random_restart or gensa")

  # Start the timer
  t0 <- Sys.time()

  # Store old random state
  if (!exists(".Random.seed", mode="numeric", envir=globalenv())) {
    sample(NA)
  }
  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv());

  # Set a new seed for reproducibility
  set.seed(1)

  # Extract the formula
  if (missing(data))
    data <- environment(formula)
  call <- match.call()
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
  if ((shift_data == TRUE) & (g2 %in% c(1, 2, 3))) {
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
    fun <- function(b) suppressWarnings(esr_rho_lp(b = b, y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
  } else if (target == "psi") {
    fun <- function(b) suppressWarnings(esr_psi_lp(b = b, y = y, x = x, alpha = alpha, g1 = g1, g2 = g2))
  }

  # Optimize the model
  if (k == 1) {
    fit <- stats::optim(par = b0, fn = fun, method = "Nelder-Mead", control = optim_ctrol)
    method <- "direct_optimization"
  } else {
    if (method == "random_restart") {
      # Evaluate the N random starting points
      noise <- matrix(stats::rnorm(2*k * random_restart_ctrl$N, sd = random_restart_ctrl$sd),
                      nrow = random_restart_ctrl$N)
      rand_start <- sweep(noise, MARGIN = 2, b0, "+")
      set.seed(Sys.time())  # Reset the seed

      # If G2 is 1, 2, 3 and we do not shift the data, ensure that x'be < 0 by shifting the
      # intercepts of the starting values
      if ((shift_data == FALSE) & (g2 %in% c(1, 2, 3))) {
        max_xe <- sapply(1:nrow(rand_start), function(i) {
          max(x %*% rand_start[i, (k+1):(2*k)])
        })
        rand_start[,k+1] <- rand_start[,k+1] - (max_xe + 0.1) * (max_xe > 0)
      }

      # Evalaute losses
      rand_start <- cbind(apply(rand_start, 1, fun), rand_start)

      # Optimize over the M best and the original starting value
      start_values <- rbind(rand_start[order(rand_start[, 1])[1:random_restart_ctrl$M], -1], b0)
      fits <- apply(start_values, 1, function(b)
        stats::optim(par = b, fn = fun, method = "Nelder-Mead", control = optim_ctrol))

      # Find the best fit
      fit <- fits[[which.min(sapply(fits, "[[", "value"))]]
    } else if (method == "gensa") {
      if (!("GenSA" %in% rownames(utils::installed.packages()))) {
        stop("GenSA needed for this function to work. Please install it.")
      }
      fit <- GenSA::GenSA(par = b0, fn = fun,
                          lower = b0 - rep(gensa_ctrl$box, 2 * k),
                          upper = b0 + rep(gensa_ctrl$box, 2 * k),
                          control = gensa_ctrl[!(names(gensa_ctrl) %in% c("box"))])
    }
  }

  # Set names of the parameters
  b <- fit$par
  names(b) <- c(paste0("bq_", 1:k - 1), paste0("be_", 1:k - 1))

  # Undo the transformation
  if ((shift_data == TRUE) & (g2 %in% c(1, 2, 3))) {
    y <- y + max_y
    b[1] <- b[1] + max_y
    b[k + 1] <- b[k + 1] + max_y
    b0[1] <- b0[1] + max_y
    b0[k + 1] <- b0[k + 1] + max_y
  }

  # Reset random seed to the old state
  assign(".Random.seed", oldSeed, envir=globalenv());

  # Return results
  structure(list(call = call, formula = formula,
                 target = target, method = method, g1 = g1, g2 = g2, shift_data = shift_data,
                 alpha = alpha, y = y, x = x,
                 coefficients = b,
                 coefficients_q = b[1:k],
                 coefficients_e = b[(k + 1):(2 * k)],
                 value = fit$value,
                 time = Sys.time() - t0),
            class = "esreg")
}

#' @export
print.esreg <- function(x, digits = 4, ...) {
  cat("Call:\n")
  cat(deparse(x$call), "\n")
  cat("\nQuantile Coefficients:\n")
  print(format(x$coefficients_q, digits = digits), quote = FALSE)
  cat("\nExpected Shortfall Coefficients:\n")
  print(format(x$coefficients_e, digits = digits), quote = FALSE)
}

#' @export
fitted.esreg <- function(object, ...) {
  cbind(object$x %*% object$coefficients_q,
        object$x %*% object$coefficients_e)
}

#' @title Covariance of the joint (VaR, ES) estimator
#' @description Estimate the variance-covariance matrix of the joint (VaR, ES) estimator either using the asymptotic formulas or using the bootstrap.
#' @param object An esreg object
#' @param sparsity Sparsity estimator
#' \itemize{
#'   \item [iid] - Piecewise linear interpolation of the distribution
#'   \item nid - Hendricks and Koenker sandwich
#' }
#' @param cond_var Conditional truncated variance estimator
#' \itemize{
#'   \item [ind] Variance over all negative residuals
#'   \item scl_N Scaling with the Normal distribution
#'   \item scl_t Scaling with the t-distribution
#'   \item scl_st Scaling with the skewed t-distribution
#' }
#' @param bandwidth_type Bofinger, Chamberlain or Hall-Sheather
#' @param bootstrap_method
#' \itemize{
#'   \item [NULL] Use the asymptotic estimator
#'   \item iid bootstrap
#'   \item stationary bootstrap (Politis & Romano, 1994)
#' }
#' @param B Number of bootstrap iterations
#' @param block_length Average block length for the stationary bootstrap
#' @param ... additional arguments
#' @export
vcov.esreg <- function(object, sparsity = "iid", cond_var = "ind", bandwidth_type = "Hall-Sheather",
                       bootstrap_method = NULL, B = 1000, block_length = NULL, ...) {
  fit <- object

  if(is.null(bootstrap_method)) {
    if(!(sparsity %in% c("iid", "nid")))
      stop("sparsity can be iid or nid")
    if(!(cond_var %in% c("ind", "scl_N", "scl_t", "scl_st")))
      stop("cond_var can be ind, scl_N, scl_t or scl_st")
    if(!(bandwidth_type %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
      stop("bandwidth_type can be Bofinger, Chamberlain or Hall-Sheather")

    # Extract some elements
    y <- fit$y
    x <- fit$x
    n <- nrow(x)
    k <- ncol(x)
    coefficients_q <- fit$coefficients_q
    coefficients_e <- fit$coefficients_e

    # Transform the data and coefficients
    if ((fit$shift_data == TRUE) & (fit$g2 %in% c(1, 2, 3))) {
      max_y <- max(y)
      y <- y - max_y
      coefficients_q[1] <- coefficients_q[1] - max_y
      coefficients_e[1] <- coefficients_e[1] - max_y
    }

    # Precompute some quantities
    xq <- as.numeric(x %*% coefficients_q)
    xe <- as.numeric(x %*% coefficients_e)
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
    cv <- conditional_truncated_variance(y = u, x = x, approach = cond_var)

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
                    alpha = fit$alpha, g1 = fit$g1, g2 = fit$g2,
                    method = 'random_restart',
                    random_restart_ctrl = list(M = 1, N = 1, sd=0))
      fitb$coefficients
    })

    # Compute the covariance
    cov <- stats::cov(t(b))
  }

  # Set the names
  rownames(cov) <- colnames(cov) <- names(stats::coef(fit))

  # Return the estimated covariance
  cov
}

#' @title esreg summary
#' @description Summarize details about the regression estimates.
#' @param object An esreg object
#' @param ... Accepts all parameters you can pass to \code{\link{vcov.esreg}}.
#' @export
summary.esreg <- function(object, ...) {
  fit <- object
  x <- fit$x
  n <- nrow(x)
  k <- ncol(x)
  cov <- vcov.esreg(fit, ...)
  se <- sqrt(diag(cov))
  tval <- stats::coef(fit) / se
  coef_mat <- cbind(Estimate     = stats::coef(fit),
                    `Std. Error` = se,
                    `t value`    = tval,
                    `Pr(>|t|)`   = 2 * stats::pt(abs(tval), n - 2*k, lower.tail = FALSE))

  structure(c(fit, list(cov = cov, coef_mat = coef_mat)), class = "summary.esreg")
}

#' @export
print.summary.esreg <- function(x, ...) {
  k <- length(x$coefficients_q)
  cat("Call:\n", paste0(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\nalpha: ", sprintf("%.3f", x$alpha), "\n")
  cat(.G_function_names(x$g1, x$g2)[1], "\n")
  cat(.G_function_names(x$g1, x$g2)[2], "\n")
  cat("Value: ", sprintf("%.9f", x$value), "\n")
  cat("\nQuantile Coefficients:\n")
  stats::printCoefmat(x$coef_mat[1:k,], signif.legend = FALSE)
  cat("\nExpected Shortfall Coefficients:\n")
  stats::printCoefmat(x$coef_mat[(k+1):(2*k),])
}
