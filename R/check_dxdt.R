#' Plot estimated percentage relative error as a function of dx
#'
#' @param ds named list of distance sampling data:
#'       $data: matrix with observations for each row with (transect ID, x, y, t) columns
#'       $transect: matrix with transects for each row with (transect ID, length)
#'       $aux: vector of variables (region width, region length, truncation distance, observer speed, transect type 0line 1point)
#'       $delta: vector of discretisation sizes in (space, time)
#'       $hazardfn: code for hazard function to use (see ?hazardfns), default is 1
#'       $move: 0 = 2d CDS model, 1 = 2d MDS model with estimated diffusion (must supply tag data), 2 = 2D MDS model with fixed diffusion (must supply fixed.sd)
#' @param move named list of movement data:
#'       $fixed.sd fixed diffusion paramter, must be supplied and only used when ds$move = 2
#'       $data list with matrix for each individual tagged with (x, y, t) observations in each row, required for ds$move = 1
#' @param par named values for (detection scale, detection shape, diffusion rate (ds$move == 1))
#' @param dx lower and upper limit of dxs to try 
#' @param by increment to increase dx by for each iteration
#' @param print FALSE by default, if TRUE then useful output is printed
#'
#' @return plot of dx against error and data frame of dxs and percentage relative error
#' @export
#'
check.dx <- function(ds,
                     move,
                     par,
                     dx = NULL,
                     by = NULL, 
                     print = FALSE) {
  save_ds <- ds 
  on.exit(ds <- save_ds)
  if (is.null(dx)) dx <- c(ds$delta[1] / 10, ds$delta[1] * 10)
  if (is.null(by)) by <- ds$delta[1] / 10 
  dxs <- seq(dx[1], dx[2], by = by)
  ndx <- length(dxs)
  res <- rep(0, ndx)
  for (i in 1:ndx) {
    ds$delta[1] <- dxs[i]
    res[i] <- mds.penc(ds, move, par, print)
  }
  best <- res[1]
  res <- 100*(res-best)/best
  plot(dxs, res, type = "b", pch = 20, lwd = 1.5, xlab = "dx", ylab = "Estimated Percentage Relative Error")
  return(data.frame(dx = dxs, err = res))
}

#' Plot estimated percentage relative error as a function of dx
#'
#' @param ds named list of distance sampling data:
#'       $data: matrix with observations for each row with (transect ID, x, y, t) columns
#'       $transect: matrix with transects for each row with (transect ID, length)
#'       $aux: vector of variables (region width, region length, truncation distance, observer speed, transect type 0line 1point)
#'       $delta: vector of discretisation sizes in (space, time)
#'       $hazardfn: code for hazard function to use (see ?hazardfns), default is 1
#'       $move: 0 = 2d CDS model, 1 = 2d MDS model with estimated diffusion (must supply tag data), 2 = 2D MDS model with fixed diffusion (must supply fixed.sd)
#' @param move named list of movement data:
#'       $fixed.sd fixed diffusion paramter, must be supplied and only used when ds$move = 2
#'       $data list with matrix for each individual tagged with (x, y, t) observations in each row, required for ds$move = 1
#' @param par named values for (detection scale, detection shape, diffusion rate (ds$move == 1))
#' @param dt lower and upper limit of dts to try 
#' @param by increment to increase dt by for each iteration
#' @param print FALSE by default, if TRUE then useful output is printed
#'
#' @return plot of dt against error and data frame of dts and percentage relative error
#' @export
#'
check.dt <- function(ds,
                     move,
                     par,
                     dt = NULL,
                     by = NULL, 
                     print = FALSE) {
  save_ds <- ds 
  on.exit(ds <- save_ds)
  if (is.null(dt)) dt <- c(ds$delta[2] / 10, ds$delta[2] * 10)
  if (is.null(by)) by <- ds$delta[1] / 10 
  dts <- seq(dt[1], dt[2], by = by)
  ndt <- length(dts)
  res <- rep(0, ndt)
  for (i in 1:ndt) {
    ds$delta[2] <- dts[i]
    res[i] <- mds.penc(ds, move, par, print)
  }
  best <- res[1]
  res <- 100*(res-best)/best
  plot(dts, res, type = "b", pch = 20, lwd = 1.5, xlab = "dt", ylab = "Estimated Percentage Relative Error")
  return(data.frame(dt = dts, penc = res))
}

