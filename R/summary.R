#' Summary of MDS fitted model
#'
#' @param object MDS fitted model
#' @param ... ignored 
#'
#' @return print summary; N and D are abundance and density estimates. Reported
#' abundance is for the region of area given by ds$aux[1:2]
#' @export
summary.mds <- function(object, ...) {
  mod <- object
  ds <- mod$ds
  move <- mod$move
  cat("Distance Sampling with movement model analysis\n")
  if (!is.null(mod$stratum)) cat("Stratum: ", mod$stratum, "\n")
  cat("Number of observations: ", nrow(ds$data), "\n")
  cat("Truncation distance: ", ds$aux[3], "\n\n")
  cat("Detection Model: ")
  if (ds$hazardfn == 1) cat("Isotropic radial hazard\n")
  if (ds$hazardfn == 2) cat("Anisotropic radial hazard\n")
  cat("Movement: ")
  if (ds$move == 0) cat("None\n")
  if (ds$move == 1) cat("Estimated diffusion\n")
  if (ds$move == 2) cat("Fixed diffusion: sd = ", move$fixed.sd, "\n")
  cat("\n")
  cat("Parameter Estimates:\n")
  print(round(mod$result, 4))
  cat("\n")
  cat("Mean detection probability: ", round(mod$penc, 4), "\n")
  cat("Loglik: ", -round(mod$fit$minimum, 2), "  AIC: ", round(mod$AIC, 2), "\n")
  invisible(mod)
}

#' Convert moveds s,d parameters to Hayes-Buckland b, sigma
#'
#' @param s scale parameter
#' @param d shape parameter 
#' @param v observer speed, default is one, if zero then point transect assumed
#' @param T duration of point transect
#'
#' @return named vector of shape (b) and scale (sigma) under Hayes-Buckland 
#' @export
s2sigmab <- function(s, d, v = 1, T = 1) {
  if (v > 0) {
    b <- d
    c <- s^d / v
    sigma <- c * 2 / beta(b/2, 1/2) 
    sigma <- sigma^(1/b)
  } else {
    b <- d 
    sigma <- s * T^(1/d)
  }
  return(c(b = b, sigma = sigma))
}
