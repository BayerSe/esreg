
#' Verify the derivatives of the specification functions
#'
#' Can be used to check whether the derivatives of G1 and G2_curly are properly implemented.
#'
#' @param g1 type of the g1 function
#' @param g2 type of the g2 function
#' @param z value at which we evalaut the derivatives. Random if NULL.
#' @param tol Tolerance of the derivative
.verify_g_function <- function(g1, g2, z = NULL, tol = 1e-08) {
    h <- 1e-05
    if (is.null(z)) 
        z <- -abs(stats::rnorm(1))
    diff_g1 <- (G1_fun(z + h, g1) - G1_fun(z - h, g1))/(2 * h) - G1_prime_fun(z, g1)
    diff_g2_1 <- (G2_curly_fun(z + h, g2) - G2_curly_fun(z - h, g2))/(2 * h) - G2_fun(z, g2)
    diff_g2_2 <- (G2_fun(z + h, g2) - G2_fun(z - h, g2))/(2 * h) - G2_prime_fun(z, g2)
    all(c(diff_g1, diff_g2_1, diff_g2_2) < tol)
}
