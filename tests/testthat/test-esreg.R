test_that("esreg formula interface works", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(100)

  fit <- esreg(y ~ x, alpha = 0.025)

  expect_s3_class(fit, "esreg")
  expect_named(fit$coefficients_q, c("(Intercept)", "x"))
  expect_named(fit$coefficients_e, c("(Intercept)", "x"))
  expect_length(fit$coefficients_q, 2)
  expect_length(fit$coefficients_e, 2)
})

test_that("esreg default interface works", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(100)
  xm <- matrix(x, ncol = 1, dimnames = list(NULL, "x"))

  fit <- esreg(xq = xm, xe = xm, y = y, alpha = 0.025)

  expect_s3_class(fit, "esreg")
  expect_length(fit$coefficients_q, 2)
  expect_length(fit$coefficients_e, 2)
})

test_that("esreg separate formula for Q and ES works", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  noise <- rnorm(100)
  y <- -x + (1 + 0.5 * x) * rnorm(100)

  fit <- esreg(y ~ x | x + noise, alpha = 0.025)

  expect_s3_class(fit, "esreg")
  expect_length(fit$coefficients_q, 2)
  expect_length(fit$coefficients_e, 3)
})

test_that("coef() returns combined coefficient vector", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(100)

  fit <- esreg(y ~ x, alpha = 0.025)
  coefs <- coef(fit)

  expect_length(coefs, 4)
  expect_named(coefs, c("(Intercept)", "x", "(Intercept)", "x"))
})

test_that("fitted() returns matrix with quantile and ES columns", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(100)

  fit <- esreg(y ~ x, alpha = 0.025)
  fv <- fitted(fit)

  expect_true(is.matrix(fv))
  expect_equal(nrow(fv), 100)
  expect_equal(ncol(fv), 2)
})

test_that("residuals() returns matrix with two columns", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(100)

  fit <- esreg(y ~ x, alpha = 0.025)
  res <- residuals(fit)

  expect_true(is.matrix(res))
  expect_equal(nrow(res), 100)
  expect_equal(ncol(res), 2)
})

test_that("predict() without newdata matches fitted()", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(100)

  fit <- esreg(y ~ x, alpha = 0.025)

  expect_equal(predict(fit), fitted(fit))
})

test_that("vcov() returns positive definite matrix", {
  set.seed(1)
  x <- rchisq(200, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(200)

  fit <- esreg(y ~ x, alpha = 0.025)
  cov <- vcov(fit, sparsity = "iid", sigma_est = "ind")

  expect_true(is.matrix(cov))
  expect_equal(nrow(cov), 4)
  expect_equal(ncol(cov), 4)
  expect_true(all(diag(cov) > 0))
})

test_that("summary.esreg returns correct structure", {
  set.seed(1)
  x <- rchisq(200, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(200)

  fit <- esreg(y ~ x, alpha = 0.025)
  s <- summary(fit)

  expect_s3_class(s, "summary.esreg")
  expect_true(is.matrix(s$coef_mat))
  expect_equal(nrow(s$coef_mat), 4)
  expect_equal(colnames(s$coef_mat),
               c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
})

test_that("esr_loss returns a scalar when return_mean = TRUE", {
  set.seed(1)
  r <- rnorm(100)
  q <- quantile(r, 0.025)
  e <- mean(r[r <= q])

  loss <- esr_loss(r = r, q = rep(q, 100), e = rep(e, 100), alpha = 0.025)

  expect_length(loss, 1)
  expect_true(is.finite(loss))
})

test_that("esr_loss returns a vector when return_mean = FALSE", {
  set.seed(1)
  r <- rnorm(100)
  q <- quantile(r, 0.025)
  e <- mean(r[r <= q])

  loss <- esr_loss(r = r, q = rep(q, 100), e = rep(e, 100),
                   alpha = 0.025, return_mean = FALSE)

  expect_length(loss, 100)
})

test_that("different g2 specifications produce valid estimates", {
  set.seed(1)
  x <- rchisq(100, df = 1)
  y <- -x + (1 + 0.5 * x) * rnorm(100)

  for (g2 in 1:5) {
    fit <- esreg(y ~ x, alpha = 0.025, g2 = g2)
    expect_s3_class(fit, "esreg")
    expect_length(coef(fit), 4)
  }
})
