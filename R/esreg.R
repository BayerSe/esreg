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
  mf <- stats::model.frame(formula = formula, data = data)
  x <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y <- stats::model.response(mf)

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
  structure(list(target = target, method = method, g1 = g1, g2 = g2,
                 alpha = alpha, y = y, x = x, b0 = b0,
                 par = b, par_q = b[1:k], par_e = b[(k + 1):(2 * k)],
                 value = fit$value, time = Sys.time() - t0),
            class = "esreg")
}



# Set some methods --------------------------------------------------------

#' @export
print.esreg <- function(x, ...) {
  cat("Target:", x$target, "\n")
  cat("Method:", x$method, "\n")
  cat("G1:    ", x$g1, "\n")
  cat("G2:    ", x$g2, "\n")
  cat("Value: ", x$value, "\n")
  cat("Time:  ", x$time, "\n\n")
  cat("Parameters:\n")
  print(x$par_q)
  print(x$par_e)
}

#' @export
coef.esreg <- function(object, ...) {
  object$par
}


#' Variance-Covariance Matrix
#'
#' @param object esreg object
#' @param boot If TRUE, bootstrap the covariance; else use the asymptotic formulas
#' @param ... other paramters passed to the low-level functions
#' @export
vcov.esreg <- function(object, boot = FALSE, ...) {
  if (boot) {
    esreg_covariance_boot(fit = object, ...)
  } else {
    esreg_covariance(fit = object, ...)
  }
}

#' @export
fitted.esreg <- function(object, ...) {
  cbind(object$x %*% object$par_q, object$x %*% object$par_e)
}
