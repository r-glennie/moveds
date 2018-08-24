#' Kolmogorov-Smirnov Goodness-of-fit tests   
#'
#' @param mod fitted mds model (returned from mds function)
#' @param delta use a different (dx, dt); otherwise uses mod$ds$delta
#' @param yquant quantile of forward distance to truncate at (for robustness to outliers)
#' @param nint number of integration points used to approximate PDF and CDF
#' @return plots of empirical and estimated CDFs for distances and list of outputs from ks.test function 
#' @export
mds.gof <- function(mod, delta = NULL, yquant = 0.95, nint = 1000) {
  # discretise transects and data
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
  range <- c(ds$aux[3], max(ds$data$y) + delta[1])
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
  Nperp <- floor(2 * range[1] / ds$delta[1])
  Nforw <- floor(range[2] / ds$delta[1])
  buf <- floor(ds$buffer / ds$delta[1])
  pdfm <- matrix(pdf, nrow = dis$numcells[2], Nforw)
  pdfm <- pdfm[buf:(Nperp+buf),]
  # line transect
  if (ds$aux[5] == 0) {
    x.pdf <- rowSums(pdfm)
    y.pdf <- colSums(pdfm)
    xrange <- seq(0, range[1], length = nint)
    xgrid <- seq(-ds$aux[3], ds$aux[3] - ds$delta[1], by = ds$delta[1])
    if (length(x.pdf) > length(xgrid)) x.pdf <- x.pdf[-1]
    # perp test
    gamdat <- data.frame(x = xgrid, y = x.pdf)
    k <- nrow(gamdat)
    xgam <- gam(y ~ s(x, k = k, fx = TRUE), data = gamdat) 
    x.sp0 <- function(r) {predict(xgam, newdata = data.frame(x = r))}
    x.sp <- function(x) {pmax(0, x.sp0(x))}
    int <- integrate(x.sp, 0, ds$aux[3])$value
    x.sp <- function(x) {pmax(0, x.sp0(x))/int}
    x.cdf <- function(x){sapply(x, FUN = function(i){integrate(x.sp, 0, i)$value})}
    x.ks <- ks.test(abs(ddata[,3]), "x.cdf")
    curve(x.cdf, 0, range[1], lwd = 1.5, lty = "dotted", xlab = "Perpendicular Distance", ylab = "CDF")
    plot(ecdf(abs(ddata[,3])), add = T, col = "grey80", lwd = 1.5)
    # forw test
    yrange <- seq(0, range[2], length = nint)
    ygrid <- seq(0, range[2] - ds$delta[1], by = ds$delta[1])
    range[2] <- quantile(ddata[,4], prob = yquant) 
    ddata <- ddata[ddata[,4]<=range[2],]
    yrange <- yrange[yrange <= range[2]]
    y.pdf <- y.pdf[ygrid <= range[2]]
    ygrid <- ygrid[ygrid <= range[2]]
    gamdat <- data.frame(x = ygrid, y = y.pdf)
    k <- nrow(gamdat)
    ygam <- gam(y ~ s(x, k = k, fx = TRUE), data = gamdat) 
    y.sp0 <- function(r) {predict(ygam, newdata = data.frame(x = r))}
    y.sp <- function(y) {pmax(0, y.sp0(y))}
    y.spy <- y.sp(yrange)
    int <- integrate(y.sp, 0, range[2])$value
    y.sp <- function(y) {pmax(0, y.sp0(y))/int}
    y.cdf <- function(x){sapply(x, FUN = function(i){integrate(y.sp, 0, i)$value})}
    ydat <- ddata[,4]
    y.ks <- ks.test(ydat, "y.cdf")
    curve(y.cdf, 0, range[2], lwd = 1.5, lty="dotted", xlab = "Forward Distance", ylab = "CDF")
    plot(ecdf(ddata[,4]), add = T, col = "grey80", lwd = 1.5)
    cat("Kolomgorov-Smirnov Test\n\n")
    cat("Perpendicular PDF: p-value = ", x.ks$p.value, " D = ", x.ks$statistic, "\n")
    cat("Forward PDF: p-value = ", y.ks$p.value, " D = ", y.ks$statistic, "\n\n")
    return(list(x.ks = x.ks, y.ks = y.ks))
  } else {
    # point transect
    xrange <- seq(0, range[1], length = nint)
    xgrid <- seq(-ds$aux[3], ds$aux[3] - ds$delta[1], by = ds$delta[1])
    yrange <- seq(0, range[2], length = nint)
    ygrid <- seq(0, range[2] - ds$delta[1], by = ds$delta[1])
    if (length(xgrid) * length(ygrid) < length(pdfm)) xgrid <- seq(-ds$aux[3], ds$aux[3], by = delta[1]) 
    x.vals <- rep(xgrid, length(ygrid))
    y.vals <- rep(ygrid, each = length(xgrid))
    r.vals <- sqrt(x.vals^2 + y.vals^2)
    gamdat <- data.frame(x = r.vals, y = r.vals * as.vector(pdfm))
    k <- length(xgrid)
    fitgam <- gam(y ~ s(x, k = k, fx=TRUE), data = gamdat) 
    r.sp0 <- function(r) {predict(fitgam, newdata = data.frame(x = r))}
    r.sp <- function(r) {pmax(0, r.sp0(r))}
    rrange <- seq(0,ds$aux[3])
    r.spy <- r.sp(rrange) 
    int <- integrate(r.sp, 0, ds$aux[3])$value 
    r.sp <- function(r) {pmax(0, r.sp0(r)) / int}
    r.cdf <- function(x){sapply(x, FUN = function(i){integrate(r.sp, 0, i)$value})}
    rdat <- sqrt(ddata[,3]^2 + ddata[,4]^2)
    r.ks <- ks.test(rdat, "r.cdf")
    curve(r.cdf, 0, ds$aux[3], lwd = 1.5, lty="dotted", xlab = "Radial Distance", ylab = "CDF")
    plot(ecdf(rdat), add = T, col = "grey80", lwd = 1.5)
    cat("Kolomgorov-Smirnov Test\n\n")
    cat("Radial PDF: p-value = ", r.ks$p.value, " D = ", r.ks$statistic, "\n")
    return(list(r.ks = r.ks))
  }
}
