#' Joint (VaR, ES) loss
#'
#'
#' @param r Vector of returns
#' @param q Vector of quantiles
#' @param e Vector of expected shortfall
#' @param alpha Quantile index
#' @param g1 1, 2, see \link{G1_fun}
#' @param g2 1, 2, 3, 4, 5, see \link{G2_curly_fun}, \link{G2_fun}
#' @param return_mean If TRUE returns the average tick loss, else the individual values
#' @references Fissler and Ziegel (2016)
#' @export
esr_loss <- function(r, q, e, alpha, g1 = 2L, g2 = 1L, return_mean = TRUE) {
    
    G1 <- function(z) sapply(z, G1_fun, type = g1)
    G2_curly <- function(z) sapply(z, G2_curly_fun, type = g2)
    G2 <- function(z) sapply(z, G2_fun, type = g2)
    
    loss <- ((r <= q) - alpha) * (G1(q) - G1(r)) + G2(e) * (e - q + (q - r) * (r <= q)/alpha) - G2_curly(e)
    
    if (return_mean) 
        mean(loss) else loss
}


#' Generalized piecewise linear loss
#'
#' Equivalent to the tick / check loss when g is the identity function.
#'
#' @param r Vector of returns
#' @param q Vector of quantiles
#' @param alpha Quantile
#' @param g Nondecreasing function
#' @param return_mean If TRUE returns the average tick loss, else the individual values
#'
#' @references Gneiting (2011)
#' @export
gpl <- function(r, q, alpha, g = function(x) x, return_mean = TRUE) {
    loss <- (alpha - (r <= q)) * (g(r) - g(q))
    if (return_mean) 
        mean(loss) else loss
}
