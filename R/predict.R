#' Predict abundance from fitted distance sampling model for transects
#'
#' @param object fitted mds model 
#' @param newdata if NULL then predicted abundance for all covered area, 
#'        otherwise, provide transects (in form of ds$transect described 
#'        for mds function) for transects covering stratum of interest 
#' @param stratum.name optional, used to label prediction object 
#' @param level confidence interval level, default is 0.95
#' @param ... ignored 
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
predict.mds <- function(object,
                newdata = NULL, 
                stratum.name = NULL, 
                level = 0.95,
                ...) {
  mds <- object
  mod <- mds$fit
  ds <- mds$ds
  move <- mds$move
  if (!is.null(newdata)) ds$transect <- newdata
  if (is.null(stratum.name) & !is.null(newdata)) stratum.name <- "Unnamed Stratum"  
  # only include observation made within stratum 
  ds$data <- ds$data[ds$data$transect %in% ds$transect[,1],]
  # initial parameters
  move.method <- ds$move
  # discretise transects and data
  dis <- Discretise(ds)
  dtrans <- dis$dtrans
  ddata <- dis$ddata
  fixed.sd <- 0
  if (move.method == 2) {
    fixed.sd <- move$fixed.sd
  }
  # get estimates
  estimate <- Working2Natural(mod$estimate, ds$hazardfn)
  # estimate covered area and abundance
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
  # VARIANCE ESTIMATION
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
  # cis
  alpha <- 1 - level
  lower <- Working2Natural(mod$estimate - qnorm(1 - alpha / 2) * sds, ds$hazardfn)
  upper <- Working2Natural(mod$estimate + qnorm(1 - alpha / 2) * sds, ds$hazardfn)
  A <- prod(ds$aux[1:2])
  D <- N / A
  var.D <- var.N / A^2
  N.ci <- N + c(-1, 1) * qnorm(1 - alpha / 2) * sqrt(var.N)
  D.ci <- N.ci / A
  # calc AIC
  aic <- 2 * mod$minimum + 2 * length(estimate)
  # prepare results
  avg.pdet <-penc / k
  res <- data.frame(Estimate = c(estimate, N, D),
                    SE = c(sds, sqrt(var.N), sqrt(var.D)),
                    LCL = c(lower, N.ci[1], D.ci[1]),
                    UCL = c(upper, N.ci[2], D.ci[2]))
  rownames(res) <- rownames(mds$result)
  res = list(result = res, cor = cov2cor(V), penc = avg.pdet, AIC = aic, fit = mod,
             ds = ds, move = move, level = level, stratum = stratum.name)
  class(res) <- c(class(res), "mds")
  return(res)
}

