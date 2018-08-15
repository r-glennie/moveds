#' Chi-squared Goodness-of-fit tests   
#'
#' @param mod fitted mds model (returned from mds function)
#' @param delta use a different (dx, dt); otherwise uses mod$ds$delta
#' @param yquant quantile of forward distance to truncate at (for robustness to outliers)
#' @return plots of empirical and estimated CDFs for distances and list of outputs from ks.test function 
#' @export
mds.gof <- function(mod, delta = NULL, yquant = 0.95) {
  # discretise transects and data
  ds <- mod$ds
  if (!is.null(delta)) ds$delta <- delta
  fixed.sd <- move$fixed.sd
  if (is.null(fixed.sd)) fixed.sd <- 0
  move <- mod$move
  dis <- Discretise(ds)
  dtrans <- dis$dtrans
  ddata <- dis$ddata
  move.method <- ds$move
  if (move.method != 0) move.method <- 1
  range <- c(ds$aux[3], max(ds$data$y) + ds$delta[1])
  # compute expected proportion of counts
  npar <- nrow(mod$result)
  wpar <- Natural2Working(mod$result[-npar,1], ds$hazardfn)
  # compute pdf
  pdf <- GetHist(wpar,
                 range,
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
  # project pdf onto x-y space
  Nperp <- floor(2 * range[1] / ds$delta[1])
  Nforw <- floor(range[2] / ds$delta[1])
  buf <- floor(ds$buffer / ds$delta[1])
  pdfm <- matrix(pdf, nrow = dis$numcells[2], Nforw)
  pdfm <- pdfm[buf:(Nperp+buf),]
  x.pdf <- rowSums(pdfm)
  x.pdf <- x.pdf / sum(x.pdf)
  y.pdf <- colSums(pdfm)
  y.pdf <- y.pdf / sum(y.pdf)
  # x 
  xgrid <- seq(-ds$aux[3], ds$aux[3], by = ds$delta[1])
  x.dens <- hist(ddata[,3], breaks = xgrid, plot = FALSE)$density
  if (length(x.pdf) > length(x.dens)) x.pdf <- x.pdf[-1]
  x.gof <- suppressWarnings(chisq.test(x.pdf, x.dens))
  # y   
  ygrid <- seq(min(c(0, ddata[,4])), range[2], by = ds$delta[1])
  range[2] <- quantile(ddata[,4], prob = yquant) 
  range[2] <- ds$delta[1] * floor(range[2] / ds$delta[1])
  ddata <- ddata[ddata[,4]<=range[2],]
  y.pdf <- y.pdf[ygrid <= range[2]]
  ygrid <- ygrid[ygrid <= range[2]]
  y.pdf <- y.pdf[-length(y.pdf)]
  y.dens <- hist(ddata[,4], breaks = ygrid, plot = FALSE)$density
  y.gof <- suppressWarnings(chisq.test(y.pdf, y.dens))
  cat("Pearson's Chi-squared test\n\n")
  cat("Perpendicular Distance: p-value = ", x.gof$p.value, "\n")
  cat("Forward Distance: p-value = ", y.gof$p.value, "\n\n")
  return(list(x.gof = x.gof, y.gof = y.gof))
}
