#' Discretise transects and map detections to spatial grid
#'
#' @param ds distance sampling data (see mds function)
#'
#' @return list of two matrices:
#'   $dtrans matrix with (transect ID, number of observations made, total time steps)
#'   $ddata contains row for each detection with (transect id, grid cell
#'     detection was recorded in, x distance, y distance, exact time of detection)
#'   $numcells contains number of grid cells in total space, x-distance, y-distance
#'   $ymax maximum distance in the y direction (usually longest line transect or radial distance for points)
#'   $T length of survey in time steps 
#' @export
#'
Discretise <- function(ds) {
  type <- ds$aux[5]
  # truncate data by truncation distance 
  if (type == 0) {
    ds.data <- ds$data[abs(ds$data$x)<=ds$aux[3],]
  } else {
    ds.data <- ds$data[sqrt(ds$data$x^2 + ds$data$y^2)<=ds$aux[3],]
  }
  # order detections by time 
  ds.data <- ds.data[order(ds.data[,4]),]
  # unpack variables 
  ds.trans <- ds$transect
  obs.speed <- ds$aux[4]
  numtrans <- nrow(ds$transect)
  dx <- ds$delta[1]
  dt <- ds$delta[2]
  w <- ds$aux[3]
  buf <- ds$buffer
  # compuate number of cells in x and y direction 
  Nx <- floor((2 * w + 2 * buf) / dx)
  if (type == 0) max.y <- max(ds.trans[,2])
  if (type == 1) max.y <- 2 * w
  Ny <- floor((max.y + 2 * buf) / dx)
  # dtrans: transect ID, number of observations, total time steps
  dtrans <- matrix(0, nrow = numtrans, ncol = 3)
  dtrans[,1] <- ds.trans[,1]
  dtrans[,2] <- sapply(dtrans[,1], FUN = function(i){sum(ds.data[,1]==i)})
  if (type == 0) dtrans[,3] <- floor(ds.trans[,2] / (dt * obs.speed)) - 1
  if (type == 1) dtrans[,3] <- floor(ds.trans[,2] / dt) - 1
  T <- max(dtrans[,3]) + 1
  # order by time
  dtrans <- dtrans[order(dtrans[,3]),]
  # ddata: columns are transect id, grid cell, x, y, t
  ddata <- matrix(0, nrow = nrow(ds.data), ncol = 5)
  ddata[,1] <- ds.data[,1]
  ddata[,3] <- ds.data[,2]
  ddata[,4] <- ds.data[,3]
  ddata[,5] <- ds.data[,4]
  # get grid cell for each observation
  ydist <- ds.data[,3] 
  if (type == 0) ydist <- ydist + ds.data[,4] * obs.speed
  if (type == 1) ydist <- ydist + w 
  xdist <- ds.data[,2] + w
  sy <- floor((ydist + buf) / dx)
  sx <- floor((xdist + buf) / dx)
  if (type == 0) {
    # round to nearest time inside discrete survey time 
    for (i in 1:nrow(ddata)) {
      tr <- ddata[i, 1]
      tr.t <- min(dtrans[dtrans[, 1] == tr, 3])
      if (floor(ddata[i, 5]/dt) > tr.t) ddata[i, 5] <- tr.t * dt
    }
  }
  sy[sy >= Ny] <- Ny - 1 
  sx[sx >= Nx] <- Nx - 1 
  ddata[,2] <- sx + Nx * sy
  # order ddata by time 
  ddata <- ddata[order(ddata[,5]),]
  return(list(dtrans = dtrans,
              ddata = ddata,
              numcells = c(Nx * Ny, Nx, Ny),
              ymax = max.y,
              T = T))
}
