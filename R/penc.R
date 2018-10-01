#' Compute detection probability for distance sampling with movement model
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
#' @param print FALSE by default, if TRUE then useful output is printed
#'
#' @return detection probability
#' @export
#'
mds.penc <- function(ds,
                move,
                par,
                print = FALSE) {
  if (print) cat("Checking input........")
  check.input(ds, move, par, print, 0.95)
  if (print) cat("done\n")

  # initial parameters
  move.method <- ds$move
  sdpar <- 0
  npar <- length(par)
  parnames <- names(par)
  if (is.null(parnames)) parnames <- c(paste0("par", 1:(npar-1)), "sd")
  # check inputs
  if (any(par < 0)) stop("One or more parmeters are negative. They should all be positive.")
  if (move.method == 1 & sdpar * sqrt(ds$delta[2]) > ds$delta[1] & ds$move == 1) {
    warning("For estimated diffusive rate, chosen discretisation may lead to unstable approximation.
            Reduce time-step or increase grid size.")
  }
  # discretise transects and data
  if (print) cat("Discretising space and time.......")
  dis <- Discretise(ds)
  dtrans <- dis$dtrans
  ddata <- dis$ddata
  if (print) cat("done\n")
  fixed.sd <- 0
  if (move.method == 2) {
    fixed.sd <- move$fixed.sd
  }
  # if asked to print, print headers
  # save start time
  start.time <- Sys.time();
  if (print) cat("Computing detection probability.......\n")
  ini.par <- Natural2Working(par, ds$hazardfn)
  penc <- GetPenc(ini.par,
                  dtrans,
                  ds$aux,
                  ds$delta,
                  dis$numcells,
                  dis$T,
                  dis$ymax,
                  ds$buffer,
                  fixed.sd,
                  ds$hazardfn,
                  move.method)
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time)
  if (print) cat("Computed in ", time.taken, attr(time.taken, "units"), "\n")
  return(penc)
}

