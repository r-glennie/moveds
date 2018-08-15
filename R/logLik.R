#' Return loglikelihood of a fitted MDS model
#'
#' @param object fitted MDS model
#' @param ... ignored
#'
#' @return loglikelihood with attribute df = number of parameters
#' @export
logLik.mds <- function(object, ...) {
   mod <- object
   llk <- -mod$fit$minimum
   attributes(llk)$df <- length(mod$fit$estimate)
   return(llk)
}
