#' Check inputs to mds function
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
#' @param start named initial values for (detection scale, detection shape, diffusion rate (ds$move == 1))
#' @param print FALSE by default, if TRUE then likelihood and parameters are printed after each evaluation
#' @param level confidence interval level, default is 0.95
#'
#'
#' @return invisible return, will stop with error message if inputs unacceptable
check.input <- function(ds, move, start, print, level) {
  ## check ds 
  pass <- TRUE
  if (!is.list(ds) |
      length(ds) != 7 |
      !all(names(ds) %in% c("data", "transect", "aux", "delta", "buffer", "hazardfn", "move"))) pass <- FALSE
  if (!pass) {
    stop("ds must be a list of 7 objects with names data, transect, aux, delta, buffer, hazardfn, move.")
  }
  # obs
  obs <- ds$data
  pass <- TRUE
  if (!is.data.frame(obs) |
      ncol(obs) != 4 | 
      any(colnames(obs) != c("transect", "x", "y", "t")) |
      !all(sapply(obs, FUN = is.numeric))) pass <- FALSE
  if (!pass) {
    stop("ds$obs must be a dataframe with four columns named (in this order): transect, x, y, t")
  }
  # transect
  trans <- ds$transect
  pass <- TRUE
  if (ncol(trans) != 2 | 
      !all(sapply(trans, FUN = is.numeric)) | 
      length(unique(trans[,1])) != nrow(trans)) pass <- FALSE
  if (!pass) {
    stop("ds$transect must have two columns: 1st is the UNIQUE transect ID and 
         2nd is the effort. Both must be numeric.")
  }
  # aux
  aux <- ds$aux
  pass <- TRUE
  if (!is.vector(aux) |
      length(aux) != 5 |
      any(aux < 0) | 
      !(aux[5] %in% c(0, 1))) pass <- FALSE
  if (!pass) {
    stop("ds$aux must be a vector of length 5, all positive numbers, with the final 
         entry being a 0 (line transect) or 1 (point transect).")
  }
  # delta 
  delta <- ds$delta 
  pass <- TRUE
  if (!is.vector(delta) | 
      length(delta) != 2 | 
      any(delta < 0) | 
      delta[1] > aux[1] | 
      delta[1] > aux[2] | 
      delta[2] > min(trans[,2]) /aux[4]) pass <- FALSE
  if (!pass) {
    stop("ds$delta must be a vector of length 2 with first element the spatial 
         discretisation being between 0 and the region size in the x,y directions; the
         second element is the temporal discretisation and should be between 0 and
         the time taken to survey the shortest transect.")
  }
  # buffer 
  buffer <- ds$buffer
  pass <- TRUE
  if (!is.numeric(buffer) | 
      buffer < 0) pass <- FALSE
  if (!pass) {
    stop("ds$buffer must be a positive number.")
  }
  # hazardfn 
  hzfn <- ds$hazardfn
  pass <- TRUE
  if (!is.numeric(hzfn) |
      !(hzfn %in% c(1, 2))) pass <- FALSE
  if (!pass) {
    stop("ds$hazardfn must be 1 or 2.")
  }
  # move 
  ds.move <- ds$move 
  pass <- TRUE
  if (!(ds.move %in% c(0, 1, 2))) pass <- FALSE
  if (!pass) {
    stop("ds$move must be 0, 1, or 2.")
  }
  
  ## check move 
  pass <- TRUE
  if (!is.list(move) | 
      (!("data" %in% names(move)) & 
      !("fixed.sd" %in% names(move)))) pass <- FALSE
  if (!pass) {
    stop("move must be a list with either an object named data or one named fixed.sd.")
  }
  
  # check movedata 
  if (!is.null(move$data)) {
    movedat <- move$data
    pass <- TRUE
    if (!is.list(movedat) | 
        any(sapply(movedat, ncol) != 3)) pass <- FALSE
    if (!pass) {
      stop("move$data must be a list of objects; each element of the list must have three
           columns (x,y,t), all numeric.")
    }
  }
  
  # check fixed.sd 
  if (!is.null(move$fixed.sd)) {
    sd <- move$fixed.sd 
    pass <- TRUE
    if (!is.numeric(sd) | sd < 0) pass <- FALSE
    if (!pass) {
      stop("move$fixed.sd must be a positive number.")
    }
  }
  
  ## check start 
  pass <- TRUE
  if (!is.vector(start) | 
      any(start < 0)) pass <- FALSE
  if (!pass) {
    stop("start must be a vector of positive numbers.")
  } 
  
  # check print 
  pass <- TRUE
  if (!is.logical(print)) pass <- FALSE
  if (!pass) {
    stop("print must be TRUE or FALSE.")
  }
 
  # check level
  pass <- TRUE
  if (!is.numeric(level) |
      level < 0 |
      level > 1) pass <- FALSE
  if (!pass) {
    stop("level must be a number between 0 or 1.")
  }
   
  invisible(list(ds, move, start, print, level))
}
