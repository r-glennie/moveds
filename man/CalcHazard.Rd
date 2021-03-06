% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CalcHazard}
\alias{CalcHazard}
\title{Calculates hazard of detection}
\usage{
CalcHazard(x, y, dt, observer_speed, parameter, type, hzfn)
}
\arguments{
\item{x}{relative x coordinate}

\item{y}{relative y coordinate}

\item{dt}{time increment}

\item{observer_speed}{speed of the observer}

\item{parameter}{vector of (detection shape, detection scale)}

\item{type}{transect type (0 = line, 1 = point)}

\item{hzfn}{hazard function code (see ?hazardfns)}
}
\value{
hazard of detection
}
\description{
Calculates hazard of detection
}
