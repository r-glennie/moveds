#' Hazard functions available in moveds package
#'
#' @description There are two available functions:
#'
#' 1. Isotropic Hayes & Buckland hazard of the form \deqn{h(r) = r/s^{-d}} for parameters
#' c,d
#'
#' 2. Anisotropic Hayes and Buckland variant \deqn{h(r) = (\frac{x^2}{sx^2} + \frac{y^2}{sy^2})^{-d}} for
#' parameters sx, sy, d. Number (1) is a special case where sx = sy.
#'
#' @export
hazardfns <- function() {NULL}
