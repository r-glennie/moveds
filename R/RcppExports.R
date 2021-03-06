# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Get observation location in 2D space
#'
#' @param  time time to return observer position
#' @param  strip_size size of strip in (x, y) dimensions
#' @param  buffer buffer size
#' @param  delta (dx, dt) vector
#' @param  transect_type 0 = line, 1 = point
#' @param  observer_speed speed of observer
#'
#' @return  (x, y) location of observer at time t
GetObserverPosition <- function(time, strip_size, buffer, delta, transect_type, observer_speed) {
    .Call('_moveds_GetObserverPosition', PACKAGE = 'moveds', time, strip_size, buffer, delta, transect_type, observer_speed)
}

#' Calculates the sparse transition rate matrix
#'
#' @param  num_cells vector with number of cells in (total space, x-direction,
#'     y-direction)
#' @param  sd vector of diffusive standard deviation for each behavioural state
#' @param  dx grid cell size in the space (c.f. delta(0))
#'
#' @return sparse transition rate matrix
CalcTrm <- function(num_cells, sd, dx) {
    .Call('_moveds_CalcTrm', PACKAGE = 'moveds', num_cells, sd, dx)
}

#' Diffuse probability distribution over space
#'
#' @description Calculate product of v with matrix exponential of a using
#' the Arnoldi process. Thereby diffusing the probability distribution
#' according to Brownian motion. Code is transcribed from Expokit package.
#' 
#' @note NOTICE
#' Permission to use, copy, modify, and distribute EXPOKIT and its
#'   supporting documentation for non-commercial purposes, is hereby
#'     granted without fee, provided that this permission message and
#'     copyright notice appear in all copies. Approval must be sought for
#'       commercial purposes as testimony of its usage in applications.
#'     
#'     Neither the Institution (University of Queensland) nor the Author
#'       make any representations about the suitability of this software for
#'         any purpose.  This software is provided ``as is'' without express or
#'        implied warranty.
#'       
#'       The work resulting from EXPOKIT has been published in ACM-Transactions 
#'         on Mathematical Software, 24(1):130-156, 1998.
#'       
#'       The bibtex record of the citation:
#'         
#'         ARTICLE{EXPOKIT,
#'                  AUTHOR  = {Sidje, R. B.},
#'                  TITLE   = {{Expokit.} {A} Software Package for
#'                    Computing Matrix Exponentials},
#'                    JOURNAL = {ACM Trans. Math. Softw.},
#'                    VOLUME  = {24},
#'                    NUMBER  = {1},
#'                    PAGES   = {130-156}
#'           YEAR    = {1998}
#'         }
#'       
#'       Certain elements of the current software may include inadequacies
#'         that may be corrected at any time, as they are discovered. The Web 
#'         always contains the latest updates.
#'       
#'       Original Author:
#'         Roger B. Sidje <rbs@maths.uq.edu.au>
#'         Department of Mathematics, University of Queensland 
#'         Brisbane, QLD-4072, Australia, (c) 1996-2006 All Rights Reserved
#'
#' @param a transition rate matrix
#' @param  v vector to be multiplied
#' @param  t time to diffuse over
#' @param  num_cells vector with number of cells in (total space, x-direction,
#'     y-direction)
#' @param  krylov_dim dimension of the approximating Krylov space
#' @param  tol tolerance in error
#'
#' @return  diffused probability distribution
Diffuse <- function(a, v, t, num_cells, krylov_dim = 30L, tol = 1e-10) {
    .Call('_moveds_Diffuse', PACKAGE = 'moveds', a, v, t, num_cells, krylov_dim, tol)
}

#' Calculates the initial distribution of animal locations.
#' Assumes uniform distribution relative to transect.
#'
#' @param  num_cells vector with number of cells in (total space, x-direction,
#'     y-direction)
#' @param  delta spatial and temporal increments (dx, dt)
#' @param  region_size size of survey region in (x,y) extents
#'
#' @return Row vector with i^th entry probability animal in i^th grid cell initially
CalcInitialDistribution <- function(num_cells, delta, region_size) {
    .Call('_moveds_CalcInitialDistribution', PACKAGE = 'moveds', num_cells, delta, region_size)
}

#' Transform working parameters (for the optimiser) to natural parameters
#' @param working_parameter working parameters 
#' @param hzfn hazard function type 
#' @return natural parameters 
Working2Natural <- function(working_parameter, hzfn = 1L) {
    .Call('_moveds_Working2Natural', PACKAGE = 'moveds', working_parameter, hzfn)
}

#' Transform natural parameters to unconstrained working parameters
#' @param parameter natural parameters 
#' @param hzfn hazard function type 
#' @return working parameters 
Natural2Working <- function(parameter, hzfn = 1L) {
    .Call('_moveds_Natural2Working', PACKAGE = 'moveds', parameter, hzfn)
}

#' Calculates hazard of detection
#'
#' @param  x relative x coordinate
#' @param  y relative y coordinate
#' @param  dt time increment
#' @param observer_speed speed of the observer
#' @param parameter vector of (detection shape, detection scale)
#' @param  type transect type (0 = line, 1 = point)
#' @param hzfn hazard function code (see ?hazardfns)
#'
#' @return  hazard of detection
CalcHazard <- function(x, y, dt, observer_speed, parameter, type, hzfn) {
    .Call('_moveds_CalcHazard', PACKAGE = 'moveds', x, y, dt, observer_speed, parameter, type, hzfn)
}

#' Computes the probability of survival for each spatial location
#'
#' @param t time step
#' @param parameter (scale, shape, diffusion) parameter
#' @param num_cells number of cells in (x, y, all) dimensions
#' @param delta (dx, dt) vector
#' @param strip_size size of strip in (x, y) dimensions
#' @param buffer buffer size
#' @param observer_speed speed of the observer
#' @param type transect type
#' @param hzfn hazard function code
#' @param nint not used
#'
#'  @return row vector of survival probabilities over space
CalcSurvivalPr <- function(t, parameter, num_cells, delta, strip_size, buffer, observer_speed, type, hzfn, nint = 4L) {
    .Call('_moveds_CalcSurvivalPr', PACKAGE = 'moveds', t, parameter, num_cells, delta, strip_size, buffer, observer_speed, type, hzfn, nint)
}

#' Thins probability distribution by the proportion detected in each grid cell
#'
#' @param t time step
#' @param pr probability distribution over finite grid
#' @param parameter (scale, shape, diffusion) parameter
#' @param num_cells number of cells in (x, y, all) dimensions
#' @param delta (dx, dt) vector
#' @param strip_size size of strip in (x, y) dimensions
#' @param buffer buffer width 
#' @param observer_speed speed of the observer
#' @param type transect type
#' @param hzfn hazard function code 
#'  
#'  @return  thinned probability distribution
Detect <- function(t, pr, parameter, num_cells, delta, strip_size, buffer, observer_speed, type, hzfn) {
    .Call('_moveds_Detect', PACKAGE = 'moveds', t, pr, parameter, num_cells, delta, strip_size, buffer, observer_speed, type, hzfn)
}

#' Compute hazard of each detection within time-step
#'
#' @param data (x, y, t) data matrix
#' @param dt time step
#' @param transdat transect data matrix
#' @param parameter (scale, shape, diffusion) parameters
#' @param observer_speed speed of observer
#' @param type 1 = point, 0 = line transect
#' @param hzfn hazard function code (see ?hazardfns)
#' @return PDF for within-timestep detection
CalcHazardDetected <- function(data, dt, transdat, parameter, observer_speed, type, hzfn) {
    .Call('_moveds_CalcHazardDetected', PACKAGE = 'moveds', data, dt, transdat, parameter, observer_speed, type, hzfn)
}

#' Calculates movement model log-likelihood
#'
#' @param  sd diffusion
#' @param  data Rcpp List where each component represent an individual path
#' and continas a matrix where each row is an observed location (x,y,t)
#'
#' @return log-likelihood
CalcMovementLogLikelihood <- function(sd, data) {
    .Call('_moveds_CalcMovementLogLikelihood', PACKAGE = 'moveds', sd, data)
}

#' Computes what grid cells are inside and outside transect
#'
#' @param  num_cells number of cells in (total, x, y) direction 
#' @param strip_size size of strip in (x,y) directions
#' @param dx grid cell size 
#' @param w for lines, half-width, for points radius 
#' @param ymax maximum forward distance for lines 
#' @param buffer distance
#' @param type =0 for lines, =1 for points 
#'
#' @return vector with 1 for each grid cell inside and 0 otherwise 
InTransect <- function(num_cells, strip_size, dx, w, ymax, buffer, type) {
    .Call('_moveds_InTransect', PACKAGE = 'moveds', num_cells, strip_size, dx, w, ymax, buffer, type)
}

#' Computes negative log-likelihood of moveDs model
#'
#' @param  working_parameter unconstrained version of parameter vector containing
#'     (detection shape, detection scale, diffusion sd)
#' @param start start value for parameters on natural scale 
#' @param  data matrix with (trans id, grid cell,t) distance sampling survey data (assumed to be ordered by transect and time)
#' @param  transdat matrix with (stripsize(1), numcells in y, totaltimestep, number of observations)
#' @param  auxiliary_data vector containing (area x extent, area y extent, strip width, transect_type)
#' @param  delta vector of (dx, dt) spacetime increments
#' @param  num_cells number of cells in (total space, x-direction, y-direction)
#' @param  T total time of survey for longest transect
#' @param  ymax maximum length of a transect
#' @param  buffer buffer distance
#' @param  movement_data field object where each component represents an individual
#'   path and contains a matrix where each row is an observed location (x,y,t)
#' @param fixed_sd if move_method = 2
#' @param hzfn hazard function code (see ?hazardfns)
#' @param  move_method 0 = 2d CDS model, 1 = 2d MDS model (movement estimated),
#'    2 = 2d MDS model (movement fixed)
#' @param  print if TRUE then print likelihood and parmeters after evaluation
#' @param con parameters are constrained to be between 1/con * start value and 
#' con * start value 
#'
#' @return  negative log-likelihood
NegativeLogLikelihood <- function(working_parameter, start, data, transdat, auxiliary_data, delta, num_cells, T, ymax, buffer, movement_data, fixed_sd = 0, hzfn = 1L, move_method = 1L, print = FALSE, con = 100) {
    .Call('_moveds_NegativeLogLikelihood', PACKAGE = 'moveds', working_parameter, start, data, transdat, auxiliary_data, delta, num_cells, T, ymax, buffer, movement_data, fixed_sd, hzfn, move_method, print, con)
}

#' Computes covered area for entire survey
#'
#' @param  working_parameter unconstrained version of parameter vector containing
#'     (detection shape, detection scale, diffusion sd)
#' @param  transdat matrix with (stripsize(1), numcells in y, totaltimestep, number of observations)
#' @param  auxiliary_data vector containing (area x extent, area y extent, strip width, transect_type)
#' @param  delta vector of (dx, dt) spacetime increments
#' @param  num_cells number of cells in (total space, x-direction, y-direction)
#' @param  T total time of survey for longest transect
#' @param  ymax maximum length of a transect
#' @param  buffer buffer distance
#' @param fixed_sd if move_method = 2
#' @param hzfn hazard function code (see ?hazardfns)
#' @param  move_method 0 = 2d CDS model, 1 = 2d MDS model (movement estimated),
#'    2 = 2d MDS model (movement fixed)
#'
#' @return  negative log-likelihood
#' @return covered area
#' unpack auxiliary data
GetPenc <- function(working_parameter, transdat, auxiliary_data, delta, num_cells, T, ymax, buffer, fixed_sd, hzfn, move_method) {
    .Call('_moveds_GetPenc', PACKAGE = 'moveds', working_parameter, transdat, auxiliary_data, delta, num_cells, T, ymax, buffer, fixed_sd, hzfn, move_method)
}

#' Computes PDF of observed detections for each (x,y) cell around the observer.
#'
#' @param  working_parameter unconstrained version of parameter vector containing
#'     (detection shape, detection scale, diffusion sd)
#' @param range to compute out to in x and y directions 
#' @param  transdat matrix with (stripsize(1), numcells in y, totaltimestep, number of observations)
#' @param  auxiliary_data vector containing (area x extent, area y extent, strip width, transect_type)
#' @param  delta vector of (dx, dt) spacetime increments
#' @param  num_cells number of cells in (total space, x-direction, y-direction)
#' @param  T total time of survey for longest transect
#' @param  ymax maximum length of a transect
#' @param  buffer buffer distance
#' @param fixed_sd if move_method = 2
#' @param hzfn hazard function code (see ?hazardfns)
#' @param  move_method 0 = 2d CDS model, 1 = 2d MDS model (movement estimated),
#'    2 = 2d MDS model (movement fixed)
#'
#' @return  matrix where (i,j) entry is cell i*dx perpendicular and j*dx forward of
#'   observer
#' unpack auxiliary data
GetHist <- function(working_parameter, range, transdat, auxiliary_data, delta, num_cells, T, ymax, buffer, fixed_sd = 0, hzfn = 1L, move_method = 1L) {
    .Call('_moveds_GetHist', PACKAGE = 'moveds', working_parameter, range, transdat, auxiliary_data, delta, num_cells, T, ymax, buffer, fixed_sd, hzfn, move_method)
}

#' Compute hazard of detection (Hayes-Buckland isotropic hazard)
#' 
#' @param x x distance 
#' @param y y distance 
#' @param dt time step
#' @param observer_speed observer speed
#' @param parameter detection parameters 
#' @param w truncation width 
#' @param type transect type 0 = point, 1 = line
#' @return hazard of detection  
NULL

#' Get Recorded forward distance once detection occurs 
#' 
#' @param x recorded x location 
#' @param y recorded y location
#' @param accu_hazard hazard accumulated up to that time 
#' @param u log random deviate 
#' @param parameter detection parameters 
#' 
#' @return recorded forward distance 
NULL

#' Get Recorded time once detection occurs 
#' 
#' @param x recorded x location 
#' @param y recorded y location
#' @param accu_hazard hazard accumulated up to that time 
#' @param u log random deviate 
#' @param parameter detection parameters 
#' 
#' @return recorded forward distance 
NULL

#' Simulate distance sampling survey
#' 
#' @param true_parameter (detection shape, scale, diffusion sd) 
#' @param N number of animals 
#' @param auxiliary_info (region x-extent ,region y-extent, survey time, dt, 
#' transect type (0 = point, 1 = line), observer_speed, number of transects, half width of transects)
#' @param dt time step 
#' @param move = 0 no mvoement, 1 Brownian motion 
#' @return Outputs csv data file 
SimulateDsData <- function(true_parameter, N, auxiliary_info, dt, move = 0L) {
    .Call('_moveds_SimulateDsData', PACKAGE = 'moveds', true_parameter, N, auxiliary_info, dt, move)
}

