#' Unload moveds package compiled code 
#'
#' @param libpath lib path 
#'
#' @return none 
#' @export
.onUnload <- function (libpath) {
  library.dynam.unload("moveds", libpath)
}
