esreg
=====

The goal of esreg is to simultaneously model the quantile and the
Expected Shortfall of a response variable given a set of covariates.

**Warning: This package is still under heavy development!**

Installation
------------

esreg is not on [CRAN](http://cran.r-project.org/) yet. You can install
it from GitHub: `devtools::install_github('BayerSe/esreg')`.

If you are using Windows, you require the
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) for compilation
of the codes.

Examples
--------

    # Load the esreg package
    library(esreg)

    # Simulate data from DGP-(2)
    set.seed(1)
    x <- rchisq(1000, df = 1)
    y <- -x + (1 + 0.1 * x) * rnorm(1000)

    # Estimate the model and the covariance
    fit <- esreg(y ~ x, alpha = 0.025)
    cov <- esreg_covariance(fit = fit, sparsity = "nid", cond_var = "scl_t")

References
----------

[A Joint Quantile and Expected Shortfall Regression
Framework](https://arxiv.org/abs/1704.02213)
