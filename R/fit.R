#' Fit distance sampling with movement model
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
#' @param print FALSE by default, if TRUE then useful output is printed
#' @param level confidence interval level, default is 0.95
#' @param ... additional arguments to be passed to nlm (the optimisation routine)
#'
#' @return Named list:
#'         $result table of estimated parameters, standard errors, and confidence intervals
#'         $cor estimated correlation matrix between parameter estimates
#'         $penc average probability of detection
#'         $AIC Akaike's information criterion score for the model
#'         $fit the fitted model object output by optim
#'         $ds ds argument 
#'         $move move argument 
#'         $level confidence level supplied 
#' @export
#'
mds <- function(ds,
                move,
                start,
                print = FALSE,
                level = 0.95,
                ...) {
  if (print) cat("Checking input........")
  check.input(ds, move, start, print, level)
  if (print) cat("done\n")

  # initial parameters
  move.method <- ds$move
  sdpar <- 0
  npar <- length(start)
  if (move.method == 1) {
    sdpar <- EstDiff(move$data)
    start[npar] <- sdpar
  }
  parnames <- names(start)
  if (is.null(parnames)) parnames <- c(paste0("par", 1:(npar-1)), "sd")
  # check inputs
  if (any(start < 0)) stop("One or more parmeters are negative. They should all be positive.")
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
  if (print) cat("Fitting model.......\n")
  if (print) cat("llk", "    ", "parameters", "\n")
  ini.par <- Natural2Working(start, ds$hazardfn)
  mod <- suppressWarnings(nlm(NegativeLogLikelihood,
             ini.par, 
                hessian = TRUE,
                start = start, 
                data = ddata,
                transdat = dtrans,
                auxiliary_data = ds$aux,
                delta = ds$delta,
                num_cells = dis$numcells,
                T = dis$T,
                ymax = dis$ymax,
                buffer = ds$buffer,
                movement_data = move$data,
                fixed_sd = fixed.sd,
                hzfn = ds$hazardfn,
                move_method = move.method,
                print = print,
                con = 100, 
                ...))
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time)
  if (print) cat("Model fitting completed in ", time.taken, attr(time.taken, "units"), "\n")
  # converged?
  if (mod$code != 1) warning("Model failed to converge with nlm code ", mod$code, ".")
  # get estimates
  estimate <- Working2Natural(mod$estimate, ds$hazardfn)
  # estimate covered area and abundance
  if (print) cat("Estimating abundance......")
  penc <- GetPenc(mod$estimate,
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
  n <- nrow(ddata)
  N <- n / penc
  if (print) cat("done\n")
  # VARIANCE ESTIMATION
  if (print) cat("Estimating variance......")
  # variance matrix
  V <- solve(mod$hessian)
  # sds
  sds <- sqrt(diag(V))
  # variance of n
  k <- nrow(ds$transect) # number of transects
  var.n <- k * var(dtrans[,2])
  # variance of penc using sandwich estimator
  grad.penc <- numDeriv::grad(GetPenc,
                    mod$estimate,
                    transdat = dtrans,
                    auxiliary_data = ds$aux,
                    delta = ds$delta,
                    num_cells = dis$numcells,
                    T = dis$T,
                    ymax = dis$ymax,
                    buffer = ds$buffer,
                    fixed_sd = fixed.sd,
                    hzfn = ds$hazardfn,
                    move_method = move.method)
  var.penc <- t(grad.penc) %*% V %*% grad.penc
  # variance of N using delta method
  var.N <- as.numeric(N^2 * (var.n / n^2 + var.penc / penc^2))
  if (print) cat("done\n")
  # cis
  if (print) cat("Computing confidence intervals......")
  alpha <- 1 - level
  lower <- Working2Natural(mod$estimate - qnorm(1 - alpha / 2) * sds, ds$hazardfn)
  upper <- Working2Natural(mod$estimate + qnorm(1 - alpha / 2) * sds, ds$hazardfn)
  A <- prod(ds$aux[1:2])
  D <- N / A
  var.D <- var.N / A^2
  N.ci <- N + c(-1, 1) * qnorm(1 - alpha / 2) * sqrt(var.N)
  D.ci <- N.ci / A
  if (print) cat("done\n")
  # calc AIC
  if (print) cat("Preparing results......")
  aic <- 2 * mod$minimum + 2 * length(estimate)
  # prepare results
  avg.pdet <-penc / k
  res <- data.frame(Estimate = c(estimate, N, D),
                    SE = c(sds, sqrt(var.N), sqrt(var.D)),
                    LCL = c(lower, N.ci[1], D.ci[1]),
                    UCL = c(upper, N.ci[2], D.ci[2]))
  rownames(res) <- c(parnames, "N", "D")
 
  res = list(result = res, cor = cov2cor(V), penc = avg.pdet, AIC = aic, fit = mod,
             ds = ds, move = move, level = level)
  class(res) <- c(class(res), "mds")
  if (print) cat("done\n")
  return(res)
}

