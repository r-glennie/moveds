// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <errno.h>
#include <RcppArmadillo.h>
//' Get observation location in 2D space
//'
//' @param  time time to return observer position
//' @param  strip_size size of strip in (x, y) dimensions
//' @param  buffer buffer size
//' @param  delta (dx, dt) vector
//' @param  transect_type 0 = line, 1 = point
//' @param  observer_speed speed of observer
//'
//' @return  (x, y) location of observer at time t
// [[Rcpp::export]]
arma::vec GetObserverPosition(const double time,
                        const arma::vec strip_size,
                        const double buffer,
                        const arma::vec delta,
                        const int transect_type,
                        const double observer_speed) {
  arma::vec pos(2);
  pos(0) = 0.5 * strip_size(0);
  if (transect_type == 1) {
    // point transect assumed to be center of grid
	  pos(1) = 0.5 * strip_size(1);
  } else {
    pos(1) = observer_speed * time + buffer;
  }
  return pos;
}
//' Calculates the sparse transition rate matrix
//'
//' @param  num_cells vector with number of cells in (total space, x-direction,
//'     y-direction)
//' @param  sd vector of diffusive standard deviation for each behavioural state
//' @param  dx grid cell size in the space (c.f. delta(0))
//'
//' @return sparse transition rate matrix
// [[Rcpp::export]]
arma::sp_mat CalcTrm(const arma::vec num_cells, const double sd, const double dx) {
  arma::sp_mat tpr = arma::zeros<arma::sp_mat>(num_cells(0), num_cells(0));
  double rate = sd * sd / (2 * dx * dx);
  int s;
  for (int i = 0; i < num_cells(1); ++i) {
    for (int j = 0; j < num_cells(2); ++j) {
      s = i + num_cells(1) * j;
      if (i < num_cells(1) - 1) {
        tpr(s, s + 1) = rate;
      }
      if (i > 0) {
        tpr(s, s - 1) = rate;
      }
      if (j < num_cells(2) - 1) {
        tpr(s, s + num_cells(1)) = rate;
      }
      if (j > 0) {
        tpr(s, s - num_cells(1)) = rate;
      }
      tpr(s, s) = -4 * rate;
    }
  }
  return tpr.t();
}
//' Diffuse probability distribution over space
//'
//' @description Calculate product of v with matrix exponential of a using
//' the Arnoldi process. Thereby diffusing the probability distribution
//' according to Brownian motion. Code is transcribed from Expokit package.
//' 
//' @note NOTICE
//' Permission to use, copy, modify, and distribute EXPOKIT and its
//'   supporting documentation for non-commercial purposes, is hereby
//'     granted without fee, provided that this permission message and
//'     copyright notice appear in all copies. Approval must be sought for
//'       commercial purposes as testimony of its usage in applications.
//'     
//'     Neither the Institution (University of Queensland) nor the Author
//'       make any representations about the suitability of this software for
//'         any purpose.  This software is provided ``as is'' without express or
//'        implied warranty.
//'       
//'       The work resulting from EXPOKIT has been published in ACM-Transactions 
//'         on Mathematical Software, 24(1):130-156, 1998.
//'       
//'       The bibtex record of the citation:
//'         
//'         ARTICLE{EXPOKIT,
//'                  AUTHOR  = {Sidje, R. B.},
//'                  TITLE   = {{Expokit.} {A} Software Package for
//'                    Computing Matrix Exponentials},
//'                    JOURNAL = {ACM Trans. Math. Softw.},
//'                    VOLUME  = {24},
//'                    NUMBER  = {1},
//'                    PAGES   = {130-156}
//'           YEAR    = {1998}
//'         }
//'       
//'       Certain elements of the current software may include inadequacies
//'         that may be corrected at any time, as they are discovered. The Web 
//'         always contains the latest updates.
//'       
//'       Original Author:
//'         Roger B. Sidje <rbs@maths.uq.edu.au>
//'         Department of Mathematics, University of Queensland 
//'         Brisbane, QLD-4072, Australia, (c) 1996-2006 All Rights Reserved
//'
//' @param a transition rate matrix
//' @param  v vector to be multiplied
//' @param  t time to diffuse over
//' @param  num_cells vector with number of cells in (total space, x-direction,
//'     y-direction)
//' @param  krylov_dim dimension of the approximating Krylov space
//' @param  tol tolerance in error
//'
//' @return  diffused probability distribution
// [[Rcpp::export]]
arma::rowvec Diffuse(const arma::sp_mat a,
                     const arma::rowvec v,
                     const double t,
                     const arma::vec num_cells,
                     const int& krylov_dim = 30,
                     const double& tol = 1e-10) {
  double m = fmin(a.n_rows, krylov_dim);
  double anorm = norm(a, "Inf");
  double mxrej = 10;
  double mx;
  double btol = 1e-7;
  double gamma = 0.9;
  double mb = m;
  int nstep = 0;
  double t_now = 0;
  double t_step;
  double delta = 1.2;
  double t_out = fabs(t);
  double s_error = 0;
  double rndoff = anorm * 1e-16;
  int k1 = 1;
  double xm = 1 / m;
  double normv = norm(v);
  double avnorm;
  double beta = normv;
  double fact = std::pow((m + 1) / std::exp(1), m + 1) * std::sqrt(2 * M_PI * (m + 1));
  double t_new = (1.0 / anorm) * std::pow((fact * tol) / (4 * beta * anorm), xm);
  double s = std::pow(10, std::floor(std::log10(t_new)) - 1);
  t_new = std::ceil(t_new / s) * s;
  double sgn = t > 0 ? 1 : -1;
  int ireject;
  double err_loc;
  double phi1;
  double phi2;
  arma::vec w = v.t();
  double hump = normv;
  arma::mat vmat = arma::zeros<arma::mat>(a.n_rows, m + 1);
  arma::mat hmat = arma::zeros<arma::mat>(m + 2, m + 2);
  arma::mat fmat;
  arma::vec p;
  while (t_now < t_out) {
    Rcpp::checkUserInterrupt();
    ++nstep;
    t_step = fmin(t_out - t_now, t_new);
    vmat.zeros();
    hmat.zeros();
    vmat.col(0) = (1 / beta) * w;
    for (int j = 0; j < m; ++j) {
      p = a * vmat.col(j);
      for (int i = 0; i <= j; ++i) {
        hmat(i, j) = dot(vmat.col(i), p);
        p -= hmat(i, j) * vmat.col(i);
      }
      s = norm(p);
      if (s < btol) {
        k1 = 0;
        mb = j;
        t_step = t_out - t_now;
        break;
      }
      hmat(j + 1, j) = s;
      vmat.col(j + 1) = (1 / s) * p;
    }
    if (k1 != 0) {
      hmat(m + 1, m) = 1;
      avnorm = norm(a * vmat.col(m));
    }
    ireject = 0;
    while (ireject <= mxrej) {
      mx = mb + k1;
      fmat = expmat(sgn * t_step * hmat.submat(0, 0, mx, mx));

      if (k1 == 0) {
        err_loc = btol;
        break;
      }
      else {
        phi1 = fabs(beta * fmat(m, 0));
        phi2 = fabs(beta * fmat(m + 1, 0) * avnorm);
        if (phi1 > 10 * phi2) {
          err_loc = phi2;
          xm = 1 / m;
        }
        else if (phi1 > phi2) {
          err_loc = (phi1 * phi2) / (phi1 - phi2);
          xm = 1 / m;
        }
        else {
          err_loc = phi1;
          xm = 1 / (m - 1);
        }
      }
      if (err_loc <= delta * t_step * tol) break;
      else {
        t_step = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
        s = std::pow(10, std::floor(std::log10(t_step)) - 1);
        t_step = std::ceil(t_step / s) * s;
        if (ireject == mxrej) {
          Rcpp::Rcout << "error: requested tolerance too high for Krylov approximation" << std::endl;
        }
        ++ireject;
      }
    }
    mx = mb + fmax(0, k1 - 1);
    w = vmat.cols(0, mx) * beta * fmat.col(0).rows(0, mx);
    beta = norm(w);
    hump = fmax(hump, beta);

    t_now = t_now + t_step;
    t_new = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
    s = std::pow(10, std::floor(std::log10(t_new) - 1));
    t_new = std::ceil(t_new / s) * s;

    err_loc = fmax(err_loc, rndoff);
    s_error += err_loc;
  }
  double err = s_error;
  hump = hump / normv;
  return w.t();
}
//' Calculates the initial distribution of animal locations.
//' Assumes uniform distribution relative to transect.
//'
//' @param  num_cells vector with number of cells in (total space, x-direction,
//'     y-direction)
//' @param  delta spatial and temporal increments (dx, dt)
//' @param  region_size size of survey region in (x,y) extents
//'
//' @return Row vector with i^th entry probability animal in i^th grid cell initially
// [[Rcpp::export]]
arma::rowvec CalcInitialDistribution(const arma::vec num_cells,
                                     const arma::vec delta,
                                     const arma::vec region_size) {
  arma::rowvec initial_phi = arma::ones<arma::rowvec>(num_cells(0));
  initial_phi *= delta(0) * delta(0);
  initial_phi /=  prod(region_size);
  return(initial_phi);
}
//' Transform working parameters (for the optimiser) to natural parameters
//' @param working_parameter working parameters 
//' @param hzfn hazard function type 
//' @return natural parameters 
// [[Rcpp::export]]
arma::vec Working2Natural(arma::vec working_parameter, int hzfn = 1) {
  arma::vec parameter = arma::exp(working_parameter);
  //parameter(1) += 2; 
  return parameter;
}
//' Transform natural parameters to unconstrained working parameters
//' @param parameter natural parameters 
//' @param hzfn hazard function type 
//' @return working parameters 
// [[Rcpp::export]]
arma::vec Natural2Working(arma::vec parameter, int hzfn = 1) {
  arma::vec working_parameter(parameter);
  //working_parameter(1) -= 2; 
  working_parameter = arma::log(working_parameter);
  return working_parameter;
}

//' Calculates hazard of detection
//'
//' @param  x relative x coordinate
//' @param  y relative y coordinate
//' @param  dt time increment
//' @param observer_speed speed of the observer
//' @param parameter vector of (detection shape, detection scale)
//' @param  type transect type (0 = line, 1 = point)
//' @param hzfn hazard function code (see ?hazardfns)
//'
//' @return  hazard of detection
// [[Rcpp::export]]
double CalcHazard(const double x,
                  const double y,
                  const double dt,
                  const double observer_speed,
                  const arma::vec parameter,
                  const int type,
                  const int hzfn) {
  double hazard = 0;
  double r0, r1, abeta, abeta2, y1;
  double s, sx, sy, d, c, k; 
  switch(hzfn) {
  case 0:
    // Hayes and Buckland isotropic h(r) = (r/s)^(-2)
    // parameter = (s, d)
    s = parameter(0); 
    d = 1; 
    c = pow(s, d); 
    r0 = x * x + y * y;
    if (type == 1) {
      hazard = dt * c / pow(r0, 0.5 * d);
    }
    else {
      // assume cannot detect behind observer
      if (y < 0) return 0;
      abeta = 0.5 * (d - 1.0);
      y1 = y - observer_speed * dt;
      if (y1 < 0) y1 = 0;
      r1 = x * x  + y1 * y1;
      if (r1 < 1e-10) return arma::datum::inf;
      if (fabs(x) < 1e-10) {
        if (fabs(d - 1) < 1e-10) {
          hazard = log(sqrt(r1)) - log(sqrt(r0));
          hazard *= c;
        } else {
          hazard = 1.0 / pow(r1, abeta) - 1.0 / pow(r0, abeta);
          hazard *= c  / (d - 1.0);
        }
      } else {
        hazard = R::pbeta(x * x / r1, abeta, 0.5, 1, 0) - R::pbeta(x * x / r0,
                          abeta, 0.5, 1, 0);
        hazard *= R::beta(abeta, 0.5) * c / (2.0 * pow(fabs(x),
                                             d - 1.0));
      }
    }
    return hazard;
    break;
  
  case 1:
    // Hayes and Buckland isotropic h(r) = (r/s)^(-d)
    // parameter = (s, d)
    s = parameter(0); 
    d = 1 + parameter(1); 
    c = pow(s, d); 
    r0 = x * x + y * y;
    if (type == 1) {
      hazard = dt * c / pow(r0, 0.5 * d);
    }
    else {
      // assume cannot detect behind observer
      if (y < 0) return 0;
      abeta = 0.5 * (d - 1.0);
      y1 = y - observer_speed * dt;
      if (y1 < 0) y1 = 0;
      r1 = x * x  + y1 * y1;
      if (r1 < 1e-10) return arma::datum::inf;
      if (fabs(x) < 1e-10) {
        if (fabs(d - 1) < 1e-10) {
          hazard = log(sqrt(r0)) - log(sqrt(r1));
          hazard *= c;
        } else {
          hazard = 1.0 / pow(r1, abeta) - 1.0 / pow(r0, abeta);
          hazard *= c / (d - 1.0);
        }
      } else {
        hazard = R::pbeta(x * x / r1, abeta, 0.5, 1, 0) - R::pbeta(x * x / r0,
                          abeta, 0.5, 1, 0);
        hazard *= R::beta(abeta, 0.5) * c / (2.0 * pow(fabs(x),
                                                 d - 1.0)); 
      }
    }
    return hazard;
    break;

  case 2:
    // Hayes and Buckland anisotropic h(r) = (x^2/sx^2 + y^2/sy^2)^(-d/2)
    // parameter = (sx, sy, d)
    sx = parameter(0);
    sy = parameter(1); 
    d = 1 + parameter(2); 
    r0 = (x * x) / (sx * sx) + (y * y) / (sy * sy);
    if (type == 1) {
      hazard = dt * pow(r0, -0.5 * d); 
    }
    else {
      // assume cannot detect behind observer
      if (y < 0) return 0;
      abeta = 0.5 * (d - 1.0);
      y1 = y - observer_speed * dt;
      if (y1 < 0) y1 = 0;
      r1 = (x * x) / (sx * sx)  + (y1 * y1) / (sy * sy);
      if (r1 < 1e-10) return arma::datum::inf;
      if (fabs(x) < 1e-10) {
        if (fabs(d - 1) < 1e-10) {
          hazard = log(sqrt(y)) - log(sqrt(y1));
          hazard *= pow(sy, d);
        } else {
          hazard = 1.0 / pow(y1, abeta) - 1.0 / pow(y, abeta);
          hazard *= pow(sy, d) / (d - 1.0);
        }
      } else {
        hazard = R::pbeta(x * x / (sx * sx * r1), abeta, 0.5, 1, 0) - R::pbeta(x * x / (sx * sx * r0),
                          abeta, 0.5, 1, 0);
        hazard *= R::beta(abeta, 0.5) * pow(sx, d - 1) * sy / (2.0 * pow(fabs(x),
                                             d - 1.0));
      }
    }
    return hazard;
    break;
    
  case 3:
    // Hayes and Buckland shape-anisotropic h(r) = (x^2+(y+k)^2)^(-d/2)
    // parameter = (s, d, k)
    s = parameter(0); 
    d = 1 + parameter(1); 
    k = parameter(2); 
    c = pow(s, d); 
    r0 = x * x / (s * s) + (y / s + k) * (y / s + k); 
    if (type == 1) {
      hazard = dt * c * (1 + k * y / sqrt(r0)) / pow(r0, 0.5 * d);
    }
    else {
      // assume cannot detect behind observer
      if (y < 0) return 0;
      abeta = 0.5 * (d - 1.0);
      y1 = y - observer_speed * dt;
      if (y1 < 0) y1 = 0;
      r1 = x * x / (s * s) + (y1 / s + k) * (y1 / s + k); 
      if (r1 < 1e-10) return arma::datum::inf;
      if (fabs(x) < 1e-10) {
        if (fabs(d - 1) < 1e-10) {
          hazard = log(sqrt(y + s * k)) - log(sqrt(y1 + s * k));
          hazard *= c;
        } else {
          hazard = 1.0 / pow(y1 + s * k, abeta) - 1.0 / pow(y + s * k, abeta);
          hazard *= c / (d - 1.0);
        }
      } else {
        hazard = R::pbeta(x * x / (s * s * r1), abeta, 0.5, 1, 0) - R::pbeta(x * x / (s * s * r0),
                          abeta, 0.5, 1, 0);
        hazard *= R::beta(abeta, 0.5) * c / (2.0 * pow(fabs(x), d - 1.0)); 
      }
    }
    return hazard;
    break;
    
  case 4:
    // Hayes and Buckland anisotropic h(r) = (x^2/sx^2 + (y/sy + k)^2)^(-d/2)
    // parameter = (sx, sy, d, k)
    sx = parameter(0);
    sy = parameter(1); 
    d = 1 + parameter(2);
    k = parameter(3); 
    r0 = (x * x) / (sx * sx) + pow(y / sy + k, 2.0); 
    if (type == 1) {
      hazard = dt * pow(r0, -0.5 * d); 
    }
    else {
      // assume cannot detect behind observer
      if (y < 0) return 0;
      abeta = 0.5 * (d - 1.0);
      y1 = y - observer_speed * dt;
      if (y1 < 0) y1 = 0;
      r1 = (x * x) / (sx * sx)  + pow(y1 / sy + k, 2.0); 
      if (r1 < 1e-10) return arma::datum::inf;
      if (fabs(x) < 1e-10) {
        if (fabs(d - 1) < 1e-10) {
          hazard = log(sqrt(y)) - log(sqrt(y1));
          hazard *= pow(sy, d);
        } else {
          hazard = 1.0 / pow(y1, abeta) - 1.0 / pow(y, abeta);
          hazard *= pow(sy, d) / (d - 1.0);
        }
      } else {
        hazard = R::pbeta(x * x / (sx * sx * r1), abeta, 0.5, 1, 0) - R::pbeta(x * x / (sx * sx * r0),
                          abeta, 0.5, 1, 0);
        hazard *= R::beta(abeta, 0.5) * pow(sx, d - 1) * sy / (2.0 * pow(fabs(x),
                                            d - 1.0));
      }
    }
    return hazard;
    break;

  default:
    Rcpp::Rcout << "error: no hazard specified." << std::endl;
  return -arma::datum::inf;
  }
}
//' Computes the probability of survival for each spatial location
//'
//' @param t time step
//' @param parameter (scale, shape, diffusion) parameter
//' @param num_cells number of cells in (x, y, all) dimensions
//' @param delta (dx, dt) vector
//' @param strip_size size of strip in (x, y) dimensions
//' @param buffer buffer size
//' @param observer_speed speed of the observer
//' @param type transect type
//' @param hzfn hazard function code
//' @param nint not used
//'
//'  @return row vector of survival probabilities over space
// [[Rcpp::export]]
arma::rowvec CalcSurvivalPr(const int t,
                      const arma::vec parameter,
                      const arma::vec num_cells,
                      const arma::vec delta,
                      const arma::vec strip_size,
                      const double buffer,
                      const double observer_speed,
                      const int type,
                      const int hzfn,
                      const int nint = 4) {

  arma::rowvec pr_survive = arma::ones<arma::rowvec>(num_cells(0));
  arma::vec observer_position(GetObserverPosition(t * delta(1), strip_size, buffer, delta, type,
                                            observer_speed));
  double x, ix;
  double y, iy;
  int s;
  int ymin = 0; 
  if (type == 0) ymin = floor(observer_position(1) / delta(0)); 
  for (int x_cell = 0; x_cell < num_cells(1); ++x_cell) {
    for (int y_cell = ymin; y_cell < num_cells(2); ++y_cell) {
      s = x_cell + num_cells(1) * y_cell;
      pr_survive(s) = 0; 
      x = x_cell * delta(0) - observer_position(0);
      y = y_cell * delta(0) - observer_position(1); 
      for (int i = 0; i < nint; ++i) {
        //for (int j = 0; j < nint; ++j) {
          ix = x + (i * delta(0)) / nint; 
          //iy = y + (j * delta(0)) / nint; 
          pr_survive(s) += CalcHazard(ix, y, delta(1), observer_speed, parameter, type, hzfn);
       // }
      }
      //pr_survive(s) /= 1.0*nint*nint; 
      pr_survive(s) /= 1.0*nint; 
      pr_survive(s) = exp(-pr_survive(s)); 
    }
  }
  return pr_survive;
}

//' Thins probability distribution by the proportion detected in each grid cell
//'
//' @param t time step
//' @param pr probability distribution over finite grid
//' @param parameter (scale, shape, diffusion) parameter
//' @param num_cells number of cells in (x, y, all) dimensions
//' @param delta (dx, dt) vector
//' @param strip_size size of strip in (x, y) dimensions
//' @param buffer buffer width 
//' @param observer_speed speed of the observer
//' @param type transect type
//' @param hzfn hazard function code 
//'  
//'  @return  thinned probability distribution
// [[Rcpp::export]]
arma::rowvec Detect(const int t,
                    const arma::rowvec pr,
                    const arma::vec parameter,
                    const arma::vec num_cells,
                    const arma::vec delta,
                    const arma::vec strip_size,
                    const double buffer,
                    const double observer_speed,
                    const int type,
                    const int hzfn) {
  arma::rowvec pr_survive = CalcSurvivalPr(t, parameter, num_cells, delta,
		  strip_size, buffer, observer_speed, type, hzfn);
  pr_survive %= pr;
  return(pr_survive);
}

//' Compute hazard of each detection within time-step
//'
//' @param data (x, y, t) data matrix
//' @param dt time step
//' @param transdat transect data matrix
//' @param parameter (scale, shape, diffusion) parameters
//' @param observer_speed speed of observer
//' @param type 1 = point, 0 = line transect
//' @param hzfn hazard function code (see ?hazardfns)
//' @return PDF for within-timestep detection
// [[Rcpp::export]]
double CalcHazardDetected(const arma::mat data,
                          double dt,
                          arma::mat transdat,
                          arma::vec parameter,
                          double observer_speed,
                          int type,
                          const int hzfn) {
  arma::vec r2, cosang;
  arma::vec t_remaining = data.col(4) - floor((data.col(4)) / dt) * dt;
  double log_hazard = 0;
  for (int i = 0; i < data.n_rows; ++i) {
    log_hazard -= CalcHazard(data(i, 2), data(i, 3) + t_remaining(i) *
		    observer_speed, t_remaining(i), observer_speed, parameter, type, hzfn);
  }

  double s, sx, sy, d, c, k; 
  switch(hzfn) {
  case 0:
    s = parameter(0); 
    d = 1; 
    c = pow(s, d); 
    r2 = sum(data.cols(2, 3) % data.cols(2, 3), 1);
    log_hazard += arma::accu(log(c) - 0.5 * d * log(r2));
    break;
  
  case 1:
    s = parameter(0); 
    d = 1 + parameter(1); 
    c = pow(s, d); 
    r2 = sum(data.cols(2, 3) % data.cols(2, 3), 1);
    log_hazard += arma::accu(log(c) - 0.5 * d * log(r2));
   break;

  case 2:
    sx = parameter(0); 
    sy = parameter(1); 
    d = 1 + parameter(2);
    log_hazard += -0.5 * d * arma::accu(log(data.col(2) % data.col(2) / (sx * sx) + data.col(3) % data.col(3) / (sy * sy))); 
    break;

  case 3:
    s = parameter(0); 
    d = 1 + parameter(1);
    k = parameter(2); 
    r2 = data.col(2) % data.col(2) / (s * s) + pow(data.col(3) / s + k, 2.0); 
    log_hazard += arma::accu(- 0.5 * d * log(r2));
    break;
    
  case 4:
    sx = parameter(0); 
    sy = parameter(1); 
    d = 1 + parameter(2);
    k = parameter(3); 
    r2 = data.col(2) % data.col(2) / (sx * sx) + pow(data.col(3) / sy + k, 2.0); 
    log_hazard += arma::accu(- 0.5 * d * log(r2));
    break;
    
  default:
    Rcpp::Rcout << "error: no hazard specified." << std::endl;
    return -arma::datum::inf;
  }
  return log_hazard;
}
//' Calculates movement model log-likelihood
//'
//' @param  sd diffusion
//' @param  data Rcpp List where each component represent an individual path
//' and continas a matrix where each row is an observed location (x,y,t)
//'
//' @return log-likelihood
// [[Rcpp::export]]
double CalcMovementLogLikelihood(const double sd, const Rcpp::List data) {
  double log_likelihood = 0;
  arma::vec xdiff, ydiff;
  int ntags = data.size();
  arma::mat tag;
  for (int i = 0; i < ntags; ++i) {
    tag = Rcpp::as<arma::mat>(data(i));
    arma::vec xdiff = arma::diff(tag.col(0));
    arma::vec ydiff = arma::diff(tag.col(1));
    arma::vec tdiff = arma::diff(tag.col(2));
    log_likelihood -= tag.n_rows * log(2 * M_PI * sd * sd);
    log_likelihood -= accu(log(tdiff) + (xdiff % xdiff + ydiff % ydiff) / (2 * sd * sd * tdiff));
  }
  return log_likelihood;
}
//' Computes what grid cells are inside and outside transect
//'
//' @param  num_cells number of cells in (total, x, y) direction 
//' @param strip_size size of strip in (x,y) directions
//' @param dx grid cell size 
//' @param w for lines, half-width, for points radius 
//' @param ymax maximum forward distance for lines 
//' @param buffer distance
//' @param type =0 for lines, =1 for points 
//'
//' @return vector with 1 for each grid cell inside and 0 otherwise 
// [[Rcpp::export]]
arma::rowvec InTransect(const arma::vec num_cells,
                        const arma::vec strip_size, 
                        const double dx, 
                        const double w,
                        const double ymax, 
                        const double buffer, 
                        const int type) {
  arma::rowvec intrans = arma::zeros<arma::rowvec>(num_cells(0)); 
  double x, y, r; 
  double top, bot, lenx, leny; 
  if (type == 0) {
    for (int i = 0; i < num_cells(1); ++i) {
      x = i * dx - strip_size(0) * 0.5;
      top = fmin(x + dx, w);
      bot = fmax(x, -w);
      lenx = top - bot; 
      if (lenx < 0) lenx = 0; 
      for (int j = 0; j < num_cells(2); ++j) {
        y = j * dx - buffer;
        top = fmin(y + dx, ymax); 
        bot = fmax(y, 0); 
        leny = top - bot; 
        if (leny < 0) leny = 0; 
        intrans(i + j * num_cells(1)) = lenx * leny / (dx * dx);
      }
    }
  } else {
    for (int i = 0; i < num_cells(1); ++i) {
      x = i * dx - strip_size(0) * 0.5; 
      for (int j = 0; j < num_cells(2); ++j) {
        y = j * dx - strip_size(1) * 0.5; 
        r = x * x + y * y; 
        if (r <= w*w) intrans(i + j * num_cells(1)) = 1; 
      }
    }
  }
  return intrans; 
}


//' Computes negative log-likelihood of moveDs model
//'
//' @param  working_parameter unconstrained version of parameter vector containing
//'     (detection shape, detection scale, diffusion sd)
//' @param start start value for parameters on natural scale 
//' @param  data matrix with (trans id, grid cell,t) distance sampling survey data (assumed to be ordered by transect and time)
//' @param  transdat matrix with (stripsize(1), numcells in y, totaltimestep, number of observations)
//' @param  auxiliary_data vector containing (area x extent, area y extent, strip width, transect_type)
//' @param  delta vector of (dx, dt) spacetime increments
//' @param  num_cells number of cells in (total space, x-direction, y-direction)
//' @param  T total time of survey for longest transect
//' @param  ymax maximum length of a transect
//' @param  buffer buffer distance
//' @param  movement_data field object where each component represents an individual
//'   path and contains a matrix where each row is an observed location (x,y,t)
//' @param fixed_sd if move_method = 2
//' @param hzfn hazard function code (see ?hazardfns)
//' @param  move_method 0 = 2d CDS model, 1 = 2d MDS model (movement estimated),
//'    2 = 2d MDS model (movement fixed)
//' @param  print if TRUE then print likelihood and parmeters after evaluation
//' @param con parameters are constrained to be between 1/con * start value and 
//' con * start value 
//'
//' @return  negative log-likelihood
// [[Rcpp::export]]
double NegativeLogLikelihood(const arma::vec working_parameter,
                             const arma::vec start, 
                             const arma::mat data,
                             const arma::mat transdat,
                             const arma::vec auxiliary_data,
                             const arma::vec delta,
                             const arma::vec num_cells,
                             const int T,
                             const double ymax,
                             const double buffer,
                             const Rcpp::List movement_data,
                             const double fixed_sd = 0,
                             const int hzfn = 1,
                             const int move_method = 1,
                             const bool print = false,
                             const double con = 100) {
  // unpack auxiliary data
  arma::vec region_size(auxiliary_data.rows(0, 1));
  arma::vec strip_size(2);
  strip_size(0) = 2 * auxiliary_data(2) + 2 * buffer;
  strip_size(1) = ymax + 2 * buffer;
  double observer_speed = auxiliary_data(3);
  int num_transects = transdat.n_rows;
  int transect_type = auxiliary_data(4);
  double dx = delta(0);
  double dt = delta(1);
  // unpack parameters
  arma::vec parameter = Working2Natural(working_parameter, hzfn);
  int npar = parameter.n_elem;
  // constraints 
  for (int p = 0; p < npar; ++p) {
    if (parameter(p) > con * start(p)) return arma::datum::inf; 
    if (parameter(p) < start(p) / con) return arma::datum::inf; 
  }
  double sd = 0;
  if (move_method == 1) sd = parameter(npar - 1);
  if (move_method == 2) sd = fixed_sd;
  // setup variables
  int curtrans = 0;
  int curobs = 0;
  double pr_survived;
  double pr_outside;
  double accu_hazard;
  double pdet; 
  double llk = 0;
  // calculate initial probability in each grid cell
  arma::rowvec pr_t = CalcInitialDistribution(num_cells, delta, region_size);
  arma::rowvec old_pr_t(pr_t);
  // probability outside buffer region at t = 0
  pr_outside = 1.0 - arma::prod(strip_size) / arma::prod(region_size);
  // compute movement matrices for survey
  arma::sp_mat trm;
  arma::rowvec flux = arma::ones<arma::rowvec>(num_cells(0));
  if (move_method > 0) {
    trm = CalcTrm(num_cells, sd, dx);
    flux = Diffuse(trm.t(), flux, dt, num_cells);
    flux = 1.0 - flux;
  }
  double num_boundary_states = floor(prod(region_size) / (dx * dx)) - num_cells(0);
  // intialise variables
  double curt = floor((data(curobs, 4)) / dt);
  if (curt < 0) curt = 0; 
  if ((curt > T - 1) & (curt < T + 1)) curt = T - 1; 
  int endtime = transdat(curtrans, 2);
  arma::rowvec intrans = InTransect(num_cells, strip_size, dx, auxiliary_data(2), ymax, buffer, transect_type); 
  accu_hazard = 0;
  pdet = 0; 
  double diff = 0; 
  // compute HMM approximation
  for (int t = 0; t < T; ++t) {
    Rcpp::checkUserInterrupt();
    // add to pr_obs, the observations that occur during time interval t
    while (t == curt) {
      if (data(curobs, 1) > num_cells(0)) {
        Rcpp::Rcout << "Warning: buffer region too small to include all detections." << std::endl;
      } else {
        llk += log(pr_t(data(curobs, 1))) + accu_hazard - log(dx * dx);
      }
      ++curobs;
      if (curobs > data.n_rows - 1) {
        curt = T + 1;
      } else {
        curt = floor((data(curobs, 4)) / dt);
        if (curt < 0) curt = 0;
        if ((curt > T - 1) & (curt < T + 1)) curt = T - 1; 
      }
    }
    // thin pr_t by those that are detected
    diff = 0; 
    diff += arma::accu(pr_t % intrans); 
    pr_t = Detect(t, pr_t, parameter, num_cells, delta, strip_size, buffer, observer_speed, transect_type, hzfn);
    diff -= arma::accu(pr_t % intrans); 
    pdet += diff * exp(accu_hazard); 
    if (arma::accu(pr_t) < 1e-10) return(arma::datum::inf);
    // move animals 
    if (move_method > 0) {
      old_pr_t = pr_t;
      //move animals that are inside strip
      try {
        pr_t = Diffuse(trm, pr_t, dt, num_cells);
      } catch(...) {
        return arma::datum::inf; 
      }
      //move animals outside strip that come into strip
      pr_t += flux * pr_outside / num_boundary_states;
      //account for transversal of boundary
      pr_outside += accu(old_pr_t % flux) - accu(flux) * pr_outside / num_boundary_states;
    }
    // add contribution to accu_hazard
    pr_survived = accu(pr_t) + pr_outside;
    accu_hazard += log(pr_survived);
    // scale to avoid underflow
    pr_t /= pr_survived;
    pr_outside /=  pr_survived;
    // if transect ends divide by conditional probability
    while (endtime == t) {
     llk -= transdat(curtrans, 1) * log(pdet);
     ++curtrans;
     if (curtrans > transdat.n_rows - 1) {
       endtime = T + 1;
     } else {
       endtime = transdat(curtrans, 2);
     }
    }
 }
 // add hazard of detections
 llk += CalcHazardDetected(data, dt, transdat, parameter, observer_speed, transect_type, hzfn);
 double movement_log_likelihood = 0;
 if (move_method == 1) movement_log_likelihood = CalcMovementLogLikelihood(sd, movement_data);
 double negative_log_likelihood = -llk - movement_log_likelihood;
 if (print) {
   int old_precision = Rcpp::Rcout.precision();
   Rcpp::Rcout.precision(4);
   Rcpp::Rcout << -negative_log_likelihood << "   ";
   for (int par = 0; par < npar; ++par) Rcpp::Rcout << parameter(par) << "   ";
   Rcpp::Rcout << std::endl;
 }
 return negative_log_likelihood;
}
//' Computes covered area for entire survey
//'
//' @param  working_parameter unconstrained version of parameter vector containing
//'     (detection shape, detection scale, diffusion sd)
//' @param  transdat matrix with (stripsize(1), numcells in y, totaltimestep, number of observations)
//' @param  auxiliary_data vector containing (area x extent, area y extent, strip width, transect_type)
//' @param  delta vector of (dx, dt) spacetime increments
//' @param  num_cells number of cells in (total space, x-direction, y-direction)
//' @param  T total time of survey for longest transect
//' @param  ymax maximum length of a transect
//' @param  buffer buffer distance
//' @param fixed_sd if move_method = 2
//' @param hzfn hazard function code (see ?hazardfns)
//' @param  move_method 0 = 2d CDS model, 1 = 2d MDS model (movement estimated),
//'    2 = 2d MDS model (movement fixed)
//'
//' @return  negative log-likelihood
//' @return covered area
//' unpack auxiliary data
// [[Rcpp::export]]
double GetPenc(const arma::vec working_parameter,
               const arma::mat transdat,
               const arma::vec auxiliary_data,
               const arma::vec delta,
               const arma::vec num_cells,
               const int T,
               const double ymax,
               const double buffer,
               const double fixed_sd,
               const int hzfn,
               int move_method) {

  arma::vec region_size(auxiliary_data.rows(0, 1));
  arma::vec strip_size(2);
  strip_size(0) = 2 * auxiliary_data(2) + 2 * buffer;
  strip_size(1) = ymax + 2 * buffer;
  double observer_speed = auxiliary_data(3);
  int num_transects = transdat.n_rows;
  int transect_type = auxiliary_data(4);
  double dx = delta(0);
  double dt = delta(1);
  // unpack parameters
  arma::vec parameter = Working2Natural(working_parameter, hzfn);
  int npar = parameter.n_elem;
  double sd = 0;
  if (move_method == 1) sd = parameter(npar - 1);
  if (move_method == 2) sd = fixed_sd;
  // setup variables
  int curtrans = 0;
  int curobs = 0;
  double pr_survived;
  double pr_outside;
  double accu_hazard;
  double pdet; 
  arma::vec penc(num_transects); penc.zeros();
  // calculate initial probability in each grid cell
  arma::rowvec pr_t = CalcInitialDistribution(num_cells, delta, region_size);
  arma::rowvec old_pr_t(pr_t);
  // probability outside buffer region at t = 0
  pr_outside = 1.0 - arma::prod(strip_size) / arma::prod(region_size);
  // compute movement matrices for survey
  arma::sp_mat trm;
  arma::rowvec flux = arma::ones<arma::rowvec>(num_cells(0));
  if (move_method > 0) {
    trm = CalcTrm(num_cells, sd, dx);
    flux = Diffuse(trm.t(), flux, dt, num_cells);
    flux = 1.0 - flux;
  }
  double num_boundary_states = floor(prod(region_size) / (dx * dx)) - num_cells(0);
  // intialise variables
  int endtime = transdat(curtrans, 2);
  arma::rowvec intrans = InTransect(num_cells, strip_size, dx, auxiliary_data(2), ymax, buffer, transect_type); 
  accu_hazard = 0;
  pdet = 0; 
  double diff; 
  // compute HMM approximation
  for (int t = 0; t < T; ++t) {
    Rcpp::checkUserInterrupt();
    // thin pr_t by those that are detected
    diff = 0; 
    diff += arma::accu(pr_t % intrans); 
    pr_t = Detect(t, pr_t, parameter, num_cells, delta, strip_size, buffer, observer_speed, transect_type, hzfn);
    diff -= arma::accu(pr_t % intrans); 
    pdet += diff * exp(accu_hazard); 
    // move animals
    if (move_method > 0) {
      old_pr_t = pr_t;
      // move animals that are inside strip
      pr_t = Diffuse(trm, pr_t, dt, num_cells);
      // move animals outside strip that come into strip
      pr_t += flux * pr_outside / num_boundary_states;
      // account for transversal of boundary
      pr_outside += accu(old_pr_t % flux) - accu(flux) * pr_outside / num_boundary_states;
    }
    // add contribution to accu_hazard
    pr_survived = accu(pr_t) + pr_outside;
    accu_hazard += log(pr_survived);
    // scale to avoid underflow
    pr_t /= pr_survived;
    pr_outside /=  pr_survived;
    // if transect ends divide by conditional probability
    while (endtime == t) {
      penc(curtrans) = pdet;
      ++curtrans;
      if (curtrans > transdat.n_rows - 1) {
        endtime = T + 1;
      } else {
        endtime = transdat(curtrans, 2);
      }
    }
  }
  return arma::accu(penc);
}
//' Computes PDF of observed detections for each (x,y) cell around the observer.
//'
//' @param  working_parameter unconstrained version of parameter vector containing
//'     (detection shape, detection scale, diffusion sd)
//' @param range to compute out to in x and y directions 
//' @param  transdat matrix with (stripsize(1), numcells in y, totaltimestep, number of observations)
//' @param  auxiliary_data vector containing (area x extent, area y extent, strip width, transect_type)
//' @param  delta vector of (dx, dt) spacetime increments
//' @param  num_cells number of cells in (total space, x-direction, y-direction)
//' @param  T total time of survey for longest transect
//' @param  ymax maximum length of a transect
//' @param  buffer buffer distance
//' @param fixed_sd if move_method = 2
//' @param hzfn hazard function code (see ?hazardfns)
//' @param  move_method 0 = 2d CDS model, 1 = 2d MDS model (movement estimated),
//'    2 = 2d MDS model (movement fixed)
//'
//' @return  matrix where (i,j) entry is cell i*dx perpendicular and j*dx forward of
//'   observer
//' unpack auxiliary data
// [[Rcpp::export]]
arma::mat GetHist(const arma::vec working_parameter,
               const arma::vec range,
               const arma::mat transdat,
               const arma::vec auxiliary_data,
               const arma::vec delta,
               const arma::vec num_cells,
               const int T,
               const double ymax,
               const double buffer,
               const double fixed_sd = 0,
               const int hzfn = 1,
               int move_method = 1) {

  arma::vec region_size(auxiliary_data.rows(0, 1));
  arma::vec strip_size(2);
  strip_size(0) = 2 * auxiliary_data(2) + 2 * buffer;
  strip_size(1) = ymax + 2 * buffer;
  double observer_speed = auxiliary_data(3);
  int num_transects = transdat.n_rows;
  int transect_type = auxiliary_data(4);
  double dx = delta(0);
  double dt = delta(1);
  // unpack parameters
  arma::vec parameter = Working2Natural(working_parameter, hzfn);
  int npar = parameter.n_elem;
  double sd = 0;
  if (move_method == 1) sd = parameter(npar - 1);
  if (move_method == 2) sd = fixed_sd;
  // setup variables
  int curtrans = 0;
  int curobs = 0;
  double pr_survived;
  double pr_outside;
  double accu_hazard;
  arma::vec obspos;
  int sobs;
  int smax;
  int Nperp = floor(2 * range(0) / dx);
  int Nforw = floor(range(1) / dx);
  arma::rowvec count(num_cells(1) * Nforw); count.zeros();
  arma::rowvec accum(num_cells(1) * Nforw); accum.zeros();
  // nalive is the number of transect still being surveyed in the meta-transect
  int nalive = num_transects;
  // calculate initial probability in each grid cell
  arma::rowvec pr_t = CalcInitialDistribution(num_cells, delta, region_size);
  arma::rowvec old_pr_t(pr_t);
  // probability outside buffer region at t = 0
  pr_outside = 1.0 - arma::prod(strip_size) / arma::prod(region_size);
  // compute movement matrices for survey
  arma::sp_mat trm;
  arma::rowvec flux = arma::ones<arma::rowvec>(num_cells(0));
  if (move_method > 0) {
    trm = CalcTrm(num_cells, sd, dx);
    flux = Diffuse(trm.t(), flux, dt, num_cells);
    flux = 1.0 - flux;
  }
  double num_boundary_states = floor(prod(region_size) / (dx * dx)) - num_cells(0);
  arma::rowvec intrans = InTransect(num_cells, strip_size, dx, auxiliary_data(2), ymax, buffer, transect_type); 
  // intialise variables
  int endtime = transdat(curtrans, 2);
  accu_hazard = 0;
  // compute HMM approximation
  for (int t = 0; t < T; ++t) {
    Rcpp::checkUserInterrupt();
    obspos = GetObserverPosition(t * dt, strip_size, buffer, delta, transect_type, observer_speed);
    // observer grid cell
    sobs = num_cells(1) * floor(obspos(1) / dx);
    // sum from that point to 2 * Nperp * Nforw state or max size
    smax = sobs + num_cells(1) * Nforw - 1;
    if (smax >= num_cells(0)) smax = num_cells(0) - 1;
    if (sobs >= num_cells(0)) sobs = num_cells(0) - 1;
    // add in all animals present in cell
    count.cols(0, smax - sobs) += exp(log(pr_t.cols(sobs, smax) % intrans.cols(sobs, smax)) + accu_hazard);
    // thin pr_t by those that are detected
    pr_t = Detect(t, pr_t, parameter, num_cells, delta, strip_size, buffer, observer_speed, transect_type, hzfn);
    // subtract those animals still present (failed to be detected)
    count.cols(0, smax - sobs) -= exp(log(pr_t.cols(sobs, smax) % intrans.cols(sobs, smax)) + accu_hazard);
    if (move_method > 0) {
      old_pr_t = pr_t;
      // move animals that are inside strip
      pr_t = Diffuse(trm, pr_t, dt, num_cells);
      // move animals outside strip that come into strip
      pr_t += flux * pr_outside / num_boundary_states;
      // account for transversal of boundary
      pr_outside += accu(old_pr_t % flux) - accu(flux) * pr_outside / num_boundary_states;
    }
    // add contribution to accu_hazard
    pr_survived = accu(pr_t) + pr_outside;
    accu_hazard += log(pr_survived);
    // scale to avoid underflow
    pr_t /= pr_survived;
    pr_outside /=  pr_survived;
    // if transect ends divide by conditional probability
    while (endtime == t) {
      accum += count; 
      --nalive;
      ++curtrans;
      if (curtrans > transdat.n_rows - 1) {
        endtime = T + 1;
      } else {
        endtime = transdat(curtrans, 2);
      }
    }
  }
  return accum;
}
