#' Plots fitted PDF from MDS model
#'
#' @param x fitted mds model (returned from mds function)
#' @param delta use a different (dx, dt); otherwise uses mod$ds$delta
#' @param yquant maximum y distance quantile to plot up to
#' @param nint discrete spacing to plot approximating spline with
#' @param extract if TRUE then predicted density for x-y distances 
#' @param ... ignored 
#' @export
plot.mds <- function(x, delta = NULL, yquant = 0.95, nint = 1000, extract = FALSE, ...) {
  mod <- x
  ds <- mod$ds
  if (!is.null(delta)) ds$delta <- delta
  delta <- ds$delta 
  move <- mod$move
  fixed.sd <- move$fixed.sd
  if (is.null(move$fixed.sd)) fixed.sd <- 0 
  dis <- Discretise(ds)
  dtrans <- dis$dtrans
  ddata <- dis$ddata
  move.method <- ds$move
  if (move.method != 0) move.method <- 1
  range <- c(ds$aux[3], max(ds$transect[,2]))
  # compute expected proportion of counts
  npar <- nrow(mod$result)
  wpar <- Natural2Working(mod$result[-((npar-1):npar),1], ds$hazardfn)
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
  Nperp <- floor(2 * range[1] / delta[1])
  Nforw <- floor(range[2] / delta[1])
  buf <- floor(ds$buffer / delta[1])
  pdfm <- matrix(pdf, nrow = dis$numcells[2], Nforw)
  pdfm <- pdfm[buf:(Nperp+buf),]
  if (extract) return(pdfm)
  # line transect
  if (ds$aux[5] == 0) {
    x.pdf <- rowSums(pdfm)
    y.pdf <- colSums(pdfm)
    xrange <- seq(0, range[1], length = nint)
    xgrid <- seq(-ds$aux[3] + 0.5 * ds$delta[1], ds$aux[3] - 0.5 * ds$delta[1], by = delta[1])
    # perp plot
    if (length(x.pdf) > length(xgrid)) x.pdf <- x.pdf[-1]
    gamdat <- data.frame(x = xgrid, y = log(x.pdf+1e-10))
    k <- nrow(gamdat)
    xgam <- gam(y ~ s(x, k = k, fx = TRUE), data = gamdat) 
    x.sp0 <- function(r) {predict(xgam, newdata = data.frame(x = r))}
    x.sp <- function(x) {exp(x.sp0(x))}
    int <- integrate(x.sp, 0, range[1])$value
    x.spy <- x.sp(xrange)
    x.spy <- x.spy / int
    dens <- hist(abs(ddata[,3]), breaks = seq(0, range[1] + delta[1], by = delta[1]), plot = FALSE)$density
    plotmax <- max(c(dens, x.spy))
    hist(abs(ddata[,3]),
         breaks = seq(0, range[1] + delta[1], by = delta[1]),
         prob = T,
         main = "Perpendicular Distance PDF",
         xlab = "Perpendicular Distance",
         ylab = "Probability Density",
         border = "grey40",
         xlim = c(0, range[1]),
         ylim = c(0, 1.1 * plotmax))
    rug(abs(ddata[,3]), quiet = TRUE)
    lines(xrange, x.spy, lwd = 1.5)
    # forw plot
    yrange <- seq(0, range[2], length = nint)
    ygrid <- seq(0, range[2] - ds$delta[1], by = delta[1])
    range[2] <- quantile(ddata[,4], prob = yquant) 
    ddata <- ddata[ddata[,4]<=range[2],]
    yrange <- yrange[yrange <= range[2]]
    y.pdf <- y.pdf[ygrid <= range[2]]
    ygrid <- ygrid[ygrid <= range[2]]
    gamdat <- data.frame(x = ygrid, y = log(y.pdf+1e-10))
    k <- nrow(gamdat)
    ygam <- gam(y ~ s(x, k = k, fx=TRUE), data = gamdat) 
    y.sp0 <- function(r) {predict(ygam, newdata = data.frame(x = r))}
    y.sp <- function(y) {exp(y.sp0(y))}
    y.spy <- y.sp(yrange)
    int <- integrate(y.sp, 0, range[2])$value
    y.spy <- y.spy / int
    dens <- hist(ddata[,4], breaks = seq(0, range[2] + delta[1], by = delta[1]), plot = FALSE)$density
    plotmax <- max(c(dens, y.spy))
    ymax <- quantile(ddata[,4], prob = yquant)
    hist(ddata[,4],
         breaks = seq(0, range[2] + delta[1], by = delta[1]),
         prob = T,
         main = "Forward Distance PDF",
         xlab = "Forward Distance",
         ylab = "Probability Density",
         border = "grey40",
         ylim = c(0, plotmax * 1.1),
         xlim = c(0, ymax))
    lines(yrange, y.spy, lwd = 1.5)
    rug(ddata[,4], quiet = TRUE)
  } else {
    # point transect
    xrange <- seq(0, range[1], length = nint)
    xgrid <- seq(-ds$aux[3], ds$aux[3] - ds$delta[1], by = delta[1])
    yrange <- seq(0, range[2], length = nint)
    ygrid <- seq(0, range[2] - ds$delta[1], by = delta[1])
    if (length(xgrid) * length(ygrid) < length(pdfm)) xgrid <- seq(-ds$aux[3], ds$aux[3], by = delta[1]) 
    x.vals <- rep(xgrid, length(ygrid))
    y.vals <- rep(ygrid, each = length(xgrid))
    r.vals <- sqrt(x.vals^2 + y.vals^2)
    gamdat <- data.frame(x = r.vals, y = log(r.vals) + as.vector(log(pdfm)))
    k <- length(xgrid)
    fitgam <- gam(y ~ s(x, k = k, fx=TRUE), data = gamdat) 
    r.sp0 <- function(r) {predict(fitgam, newdata = data.frame(x = r))}
    r.sp <- function(r) {exp(r.sp0(r))}
    rrange <- seq(0,ds$aux[3])
    r.spy <- r.sp(rrange) 
    int <- integrate(r.sp, 0, ds$aux[3])$value 
    r.spy <- r.spy / int
    rdat <- sqrt(ddata[,3]^2 + ddata[,4]^2)
    dens <- hist(rdat, breaks = seq(0, ds$aux[3], by = delta[1]), plot = FALSE)$density
    plotmax <- max(c(dens, r.spy))
    hist(rdat,
         breaks = seq(0, ds$aux[3], by = delta[1]),
         prob = T,
         main = "Radial Distance PDF",
         xlab = "Radial Distance",
         ylab = "Probability Density",
         border = "grey40",
         ylim = c(0, plotmax * 1.1))
    lines(rrange, r.spy, lwd = 1.5)
    rug(rdat, quiet = TRUE)
  }
  invisible(list(mod, ds, move))
}
