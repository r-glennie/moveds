% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CalcMovementLogLikelihood}
\alias{CalcMovementLogLikelihood}
\title{Calculates movement model log-likelihood}
\usage{
CalcMovementLogLikelihood(sd, data)
}
\arguments{
\item{sd}{diffusion}

\item{data}{Rcpp List where each component represent an individual path
and continas a matrix where each row is an observed location (x,y,t)}
}
\value{
log-likelihood
}
\description{
Calculates movement model log-likelihood
}
