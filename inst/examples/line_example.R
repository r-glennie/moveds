## moveds example
library(moveds)

# Simulate line transect survey -------------------------------------------

# c, d, diffusion rate
true.par <- c(5, 3, 2.5)
hzfn <- 1
# auxiliary information: region width, region length, half width, observer_speed, transect type
aux <- c(1000, 1000, 30, 1.0, 0)
# simulation time step
sim.dt <- 1
# true abundance
true.N <- 100
# simulate movement or not?
sim.move <- 1
# transects
n.transects <- 100
transdat <- matrix(0, nr = n.transects, nc = 2)
transdat[,1] <- 1:n.transects
transdat[,2] <- rep(1000, n.transects)
# Run simulation function to output data file
hb.par <- s2sigmab(true.par[1], true.par[2])
SimulateDsData(c(rev(hb.par), true.par[3]),
               true.N,
               c(1000, 1000, 60, 1000, 1000, 1, n.transects, 0),
               sim.dt,
               sim.move)



# Simulate movement data  -------------------------------------------------
num.tagged <- 5
observation.times <- seq(0, 1000, 10)
movedat <- SimulateMovementData(num.tagged, observation.times, true.par[3])

# 1D CDS model ------------------------------------------------
library(Distance)
# read in data and get into format for Distance package
obs <- read.csv("simulated_dsdata.csv", h = FALSE)
names(obs) <- c("transect", "x", "y", "t")
region.table <- data.frame(Region.Label = 1, Area = prod(aux[1:2]))
sample.table <- data.frame(Sample.Label = transdat[,1], Region.Label = 1,
                           Effort = transdat[, 2])
obs.table <- data.frame(object = 1:nrow(obs), Region.Label = 1, Sample.Label = obs[, 1])
distances <- abs(obs$x)

# fit 1D CDS model
cds1d <- ds(distances,
            truncation = aux[3],
            key = "hr",
            adjustment = NULL,
            region.table = region.table,
            sample.table = sample.table,
            obs.table = obs.table)

# get estimated detection parameters
det.est <- as.numeric(exp(cds1d$ddf$ds$par))
# plot estimated detection function
plot(cds1d)
# get estimated density
N.cds1d <- as.numeric(cds1d$dht$individuals$N$Estimate)
cds1d.gof <- ds.gof(cds1d)

# 2D CDS model ------------------------------------------------------------
# format for mds function
ds <- list(data = obs,
           transect = transdat,
           aux = aux,
           delta = c(1, 1000),
           buffer = 0,
           hazardfn = 1,
           move = 0)
move <- list(data = movedat,
             fixed.sd = 0.5)

cds2d <- mds(ds, move, start = c(s = 5, d = 3), print = TRUE)
summary(cds2d)
plot(cds2d, delta = c(5, 1))
cds2d.gof <- mds.gof(cds2d, delta = c(5, 1))
start <- Working2Natural(cds2d$fit$estimate, ds$hazardfn)

dxs <- seq(1, 0.1, -0.1)
cds2d_llks <- rep(0, length(dxs))
for (i in seq(dxs)){
ds <- list(data = obs,
           transect = transdat,
           aux = aux,
           delta = c(dxs[i], 1000),
           buffer = 0,
           hazardfn = 1,
           move = 0)
move <- list(data = movedat,
             fixed.sd = 0.5)

# fit 2D CDS model
cds2d_llk <- mds.penc(ds, move, par = c(s = 5, d = 3))
cds2d_llks[i] <- cds2d_llk
cat("dx = ", dxs[i], " llk = ", cds2d_llk, "\n")
}

plot(dxs, cds2d_llks, type = "b")


#plot(cds2d, delta = c(cds2d$ds$delta[1], 1))
#mds.gof(cds2d, delta = c(2.5, 1))

# fit using LT2d ----------------------------------------------------------
# library(LT2D)
# xd <- abs(obs$x)
# yd <- obs$y
# wd <- aux[3]
# ymax <- 100
# h2 <- function(y,x,b) {
#   par <- exp(b)
#   return((x*x/par[1]^2 + y*y/par[2]^2)^(-par[3]/2))
# }
# 
# fName <- "h1"
# dfit <- fityx(yd[xd<=wd],
#                 xd[xd<=wd],
#                 b = log(c(0.8,3)),
#                 hr = "h1",
#                 ystart = ymax,
#                 pi.x = "pi.const",
#                 logphi = NULL,
#                 w = wd,
#                 hessian = TRUE,
#                 control = list(trace = 5))
# est <- h1.to.HB(dfit$b)
# est[2] <- est[2] - 1
# est[1] <- est[1]^(1/est[2])
# n=length(dfit$dat$x)
# L=1000*n.transects
# pdet <- phat(dfit)
# Dhat=1000^2 * (n/pdet)/(2*wd*L)


# 2D MDS model ------------------------------------------------------------
ds$move <- 1
ds$delta <- c(5, 10)
ds$buffer <- 5
mds2d <- mds(ds, move, start = c(s = 5, d = 3, sd = 2.5), print = TRUE)
summary(mds2d)
#plot(mds2d)
#mds2d.gof <- mds.gof(mds2d)
#start <- Working2Natural(mds2d$fit$estimate, ds$hazardfn)
dxs <- seq(1, 0.1, -0.1)
mds2d_pencs <- rep(0, length(dxs))
for (i in seq(dxs)){
  ds <- list(data = obs,
             transect = transdat,
             aux = aux,
             delta = c(dxs[i], 10),
             buffer = 0,
             hazardfn = 1,
             move = 1)
  move <- list(data = movedat,
               fixed.sd = 0.5)
  
  # fit 2D CDS model
  mds2d_penc <- mds.penc(ds, move, par = c(s = 2.86, d = 2.41))
  mds2d_pencs[i] <- mds2d_penc
  cat("dx = ", dxs[i], " llk = ", mds2d_penc, "\n")
}


dxs <- seq(1, 0.1, -0.25)
mmods <- vector(mode = "list", length = length(dxs))
for (i in seq(dxs)) {
ds$move <- 1
ds$delta <- c(dxs[i], 1)
ds$buffer <- 5
mds2d <- mds(ds, move, start = c(s = start[1], d = start[2], sd = start[3]), print = TRUE)
mmods[[i]] <- mds2d
cat("dx = ", dxs[i], " N = ", mds2d$result[4,1], "\n")
}
summary(mds2d)
Ns <- sapply(mmods, FUN=function(x){x$result[4,1]})
llks <- sapply(mmods, FUN = logLik)
plot(dxs,Ns,type="b")
plot(dxs, llks, type = "b")
check.convergence(mds2d, print = TRUE)
plot(mds2d, delta = c(5, 1))

# ds$move <- 2 
# move$fixed.sd <- 2.5 
# mds.fix <- mds(ds, move, start = c(s = 5, d = 3), print = TRUE)
