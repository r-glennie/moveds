## moveds example
library(moveds)

# Simulate point transect survey -------------------------------------------

# c, d, diffusion rate
true.par <- c(10, 3, 2)
hzfn <- 1
# auxiliary information: region width, region length, half width, observer_speed, transect type
aux <- c(1000, 1000, 150, 0, 1)
# simulation time step
sim.dt <- 1
# true abundance
true.N <- 200
# simulate movement or not?
sim.move <- 1
# transects
n.transects <- 100
transdat <- matrix(0, nr = n.transects, nc = 2)
transdat[,1] <- 1:n.transects
transdat[,2] <- rep(3*60, n.transects) #runif(n.transects, 100, 1000)
# Run simulation function to output data file
SimulateDsData(true.par,
               true.N,
               c(1000, 1000, 300, 3*60, 3*60, 0, n.transects, 1),
               sim.dt,
               sim.move)

# Simulate movement data  -------------------------------------------------
num.tagged <- 10
observation.times <- seq(0, 1000, 1)
movedat <- SimulateMovementData(num.tagged, observation.times, true.par[3])

# 1D CDS model ------------------------------------------------
library(Distance)
# read in data and get into format for Distance package
obs <- read.csv("simulated_dsdata.csv", h = FALSE)
names(obs) <- c("transect", "x", "y", "t")
region.table <- data.frame(Region.Label = 1, Area = prod(aux[1:2]))
sample.table <- data.frame(Sample.Label = 1:n.transects, Region.Label = 1,
                           Effort = 1)
obs.table <- data.frame(object = 1:nrow(obs), Region.Label = 1, Sample.Label = obs[, 1])
distances <- sqrt(obs$x^2 + obs$y^2)

# fit 1D CDS model
cds1d <- ds(distances,
            truncation = aux[3],
            transect = "point",
            key = "hr",
            adjustment = NULL,
            region.table = region.table,
            sample.table = sample.table,
            obs.table = obs.table)

# get estimated detection parameters
det.est <- as.numeric(exp(cds1d$ddf$ds$par))
det.est[2] <- det.est[2] * transdat[1,2]^(-1/det.est[1])
# plot estimated detection function
plot(cds1d, pdf = TRUE)
# get estimated density
N.cds1d <- as.numeric(cds1d$dht$individuals$N$Estimate)
cds1d.gof <- ds.gof(cds1d)

# 2D CDS model ------------------------------------------------------------
# format for mds function
ds <- list(data = obs,
           transect = transdat,
           aux = aux,
           delta = c(50, 3*60),
           buffer = 0,
           hazardfn = 1,
           move = 0)
move <- list(data = movedat,
             fixed.sd = 0.5)

# fit 2D CDS model
cds2d <- mds(ds, move, start = c(s = 0.5, d = 3), print = TRUE)
summary(cds2d)
plot(cds2d)
cds2d.gof <- mds.gof(cds2d)

dxs <- seq(50, 5, -1)
cds2d_pencs <- rep(0, length(dxs))
for (i in seq(dxs)){
  ds$delta <- c(dxs[i], 3*60)
  # fit 2D CDS model
  cds2d_penc <- mds.penc(ds, move, par = c(s = 7.1, d = 2.8))
  cds2d_pencs[i] <- cds2d_penc
  cat("dx = ", dxs[i], " llk = ", cds2d_penc, "\n")
}
nest <- (cds2d_pencs - cds2d_penc) / cds2d_penc
plot(dxs, abs(nest), type = "b")
abline(h = 0.01, col = "red", type = "dotted")


# 2D MDS model ------------------------------------------------------------
ds$move <- 1
ds$delta <- c(5, 1)
ds$buffer <- 5
mds2d <- mds(ds, move, start = c(s = 5, d = 3, sd = 2.5), print = TRUE)
summary(mds2d)
plot(mds2d)
mds2d.gof <- mds.gof(mds2d)

dxs <- seq(5, 0.1, -0.1)
mds2d_pencs <- rep(0, length(dxs))
for (i in seq(dxs)){
  ds$delta <- c(dxs[i], 3*60)
  # fit 2D CDS model
  mds2d_penc <- mds.penc(ds, move, par = c(s = 10.5, d = 3.1))
  mds2d_pencs[i] <- mds2d_penc
  cat("dx = ", dxs[i], " llk = ", mds2d_penc, "\n")
}
nest <- (mds2d_pencs - mds2d_penc) / mds2d_penc
plot(dxs, abs(nest)*100, type = "b")
abline(h = 0.01, col = "red", type = "dotted")

# ds$move <- 2 
# move$fixed.sd <- 2.5 
# mds.fix <- mds(ds, move, start = c(s = 5, d = 3), print = TRUE)
