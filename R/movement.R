#' Simulate movement paths from Brownian motion
#'
#' @param num.animals number of tagged animals to simulate
#' @param observation.times vector of times animal locations are recorded
#' @param diffusion diffusion rate
#'
#' @return list of matrices, each matrix corresponds to an individual and has
#'         columns (xlocation, ylocation, observation time)
#' @export
#'
SimulateMovementData  <- function(num.animals, observation.times, diffusion) {
  obs  <- vector("list", num.animals)
  for (animal in seq(num.animals)) {
    n.obs  <- length(observation.times)
    obs.x  <- obs.y  <- numeric(n.obs)
    x  <- y  <- 0
    cur.obs  <- 1
    cur.t <- 0
    next.t  <- observation.times[cur.obs]
    while (cur.obs <= length(observation.times)) {
      x <- x + diffusion * sqrt(next.t - cur.t) * rnorm(1)
      y <- y + diffusion * sqrt(next.t - cur.t) * rnorm(1)
      obs.x[cur.obs] <- x
      obs.y[cur.obs] <- y
      cur.obs  <- cur.obs + 1
      cur.t <- next.t
      next.t <- observation.times[cur.obs]
    }
    obs.x  <- obs.x - min(obs.x)
    obs.y  <- obs.y - min(obs.y)
    obs.mat  <- as.matrix(data.frame(x = obs.x, y = obs.y, observation.times))
    obs[[animal]]  <- obs.mat
  }
  return(obs)
}

#' Estimate diffusion rate from tag data
#'
#' @param tags list of matrices, one for each tagged animal, with (x, y, t) rows
#'
#' @return estimated diffusion rate
#' @export
#'
EstDiff <- function(tags) {
  ini.sd <- sd(diff(tags[[1]][,1]) / sqrt(diff(tags[[1]][,3])))
  movemod <- optimize(CalcMovementLogLikelihood,
                      interval = c(0, 100 * ini.sd),
                      data = tags,
                      maximum = TRUE)
  return(movemod$maximum)
}


