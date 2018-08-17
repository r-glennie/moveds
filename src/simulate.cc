// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <errno.h>
#include <RcppArmadillo.h>

//' Compute hazard of detection (Hayes-Buckland isotropic hazard)
//' 
//' @param x x distance 
//' @param y y distance 
//' @param dt time step
//' @param observer_speed observer speed
//' @param parameter detection parameters 
//' @param w truncation width 
//' @param type transect type 0 = point, 1 = line
//' @return hazard of detection  
double CalcHazard2(double x, double y, double dt, double observer_speed, 
                   arma::vec parameter, double w = 0, int type = 0) {
  double hazard = 0; 
  if (type == 1) {
    double r = x * x + y * y; 
    hazard = dt * parameter(0) / pow(r, 0.5 * parameter(1));  
  }
  if (type == 0) {
    if (y < 0) return 0; 
    double abeta = 0.5 * (parameter(1) - 1.0); 
    double r0 = x * x  + y * y;
    double y1 = y - observer_speed * dt; 
    if (y1 < 0) y1 = 0;  
    double r1 = x * x  + y1 * y1;
    if (r1 < 1e-10) return arma::datum::inf; 
    if (fabs(x) < 1e-10) {
      if (fabs(parameter(1) - 1) < 1e-10) {
        hazard = log(sqrt(r1)) - log(sqrt(r0));
        hazard *= parameter(0); 
      } else {
        hazard = 1.0 / pow(r1, abeta) - 1.0 / pow(r0, abeta); 
        hazard *= parameter(0) / (parameter(1) - 1.0); 
      }
    } else {
      hazard = R::pbeta(x * x / r1, abeta, 0.5, 1, 0) - R::pbeta(x * x / r0, 
                        abeta, 0.5, 1, 0); 
      hazard *= R::beta(abeta, 0.5) * parameter(0) / (2.0 * pow(fabs(x), 
                                                parameter(1) - 1.0)); 
    } 
  }
  return hazard;   
}

//' Get Recorded forward distance once detection occurs 
//' 
//' @param x recorded x location 
//' @param y recorded y location
//' @param accu_hazard hazard accumulated up to that time 
//' @param u log random deviate 
//' @param parameter detection parameters 
//' 
//' @return recorded forward distance 
double GetRecordedYPosition(double x, double y, double accu_hazard, 
                            double u, arma::vec parameter) {
  double r = x * x + y * y; 
  double y1;
  double abeta =  0.5 * (parameter(1) - 1.0); 
  if (fabs(x) < 1e-10) {
    if (fabs(parameter(1) - 1) < 1e-10) {
      y1 = sqrt(exp(2.0 * (u - accu_hazard) / parameter(0) + log(r)) - x * x);       
    } else {
      y1 = sqrt(pow((u - accu_hazard) * (parameter(1) - 1.0) / parameter(0) + 
        1 / pow(r, abeta), -1.0 / abeta) - x * x); 
    }
  } else {
    y1 = R::pbeta(x * x / r, abeta, 0.5, 1, 0) + 2 * pow(fabs(x), 
                  parameter(1) - 1.0) * (u - accu_hazard) / (parameter(0) * R::beta(abeta, 0.5)); 
    y1 = R::qbeta(y1, abeta, 0.5, 1, 0); 
    y1 = fabs(x) * sqrt(1.0 / y1 - 1.0);  
  }
  return y1; 
}

//' Get Recorded time once detection occurs 
//' 
//' @param x recorded x location 
//' @param y recorded y location
//' @param accu_hazard hazard accumulated up to that time 
//' @param u log random deviate 
//' @param parameter detection parameters 
//' 
//' @return recorded forward distance 
double GetRecordedT(double x, double y, double accu_hazard, 
                    double u, arma::vec parameter) {
  double r = x * x + y * y; 
  double t_add = (u - accu_hazard) * pow(r, 0.5 * parameter(1)) / parameter(0);
  return t_add; 
}

//' Simulate distance sampling survey
//' 
//' @param true_parameter (detection shape, scale, diffusion sd) 
//' @param N number of animals 
//' @param auxiliary_info (region x-extent ,region y-extent, survey time, dt, 
//' transect type (0 = point, 1 = line), observer_speed, number of transects, half width of transects)
//' @param dt time step 
//' @param move = 0 no mvoement, 1 Brownian motion 
//' @return Outputs csv data file 
// [[Rcpp::export]]
int SimulateDsData(arma::vec true_parameter, 
                   int N, 
                   arma::vec auxiliary_info, 
                   double dt,
                   int move = 0) {
  arma::vec area = auxiliary_info.rows(0, 1); 
  double half_width = auxiliary_info(2) * 0.5; 
  double survey_time = auxiliary_info(4);
  double transect_type = auxiliary_info(7); 
  double observer_speed = auxiliary_info(5); 
  int num_transects = auxiliary_info(6); 
  arma::vec parameter(true_parameter); 
  if (transect_type == 0) {
    parameter(0) = 2.0 * pow(parameter(0), parameter(1)) / R::beta(0.5 * parameter(1), 0.5);
  } else {
    parameter(0) = pow(parameter(0), parameter(1));
  }
  if (transect_type == 0) ++parameter(1); 
  arma::vec observer_position(2);
  observer_position(0) = area(0) * 0.5; 
  observer_position(1) = area(1) * 0.5; 
  arma::vec detection_parameter(parameter.rows(0,1)); 
  double sd = parameter(2);
  // open data file 
  std::ofstream data("./simulated_dsdata.csv"); 
  arma::vec x(N);
  arma::vec y(N);  
  arma::vec u(N); 
  arma::vec detected(N); 
  arma::vec accu_hazard(N);
  double recorded_y; 
  double recorded_t; 
  double recorded_r; 
  double hazard;  
  bool include; 
  for (int transect = 0; transect < num_transects; ++transect) {
    if (transect_type == 0) observer_position(1) = 0.0;
    x = arma::randu(N) * area(0); 
    y = arma::randu(N) * area(1); 
    u = -log(arma::randu(N)); 
    detected.zeros(); 
    accu_hazard.zeros(); 
    for (int t = 0; t < floor(survey_time / dt); ++t) {
      for (int i = 0; i < N; ++i) {
        hazard = CalcHazard2(x(i) - observer_position(0), y(i) - 
          observer_position(1), dt, observer_speed, parameter, half_width, transect_type);
        if ((accu_hazard(i) + hazard >= u(i)) & (accu_hazard(i) < u(i))) {
          include = false; 
          if (transect_type == 0) {
            recorded_y = GetRecordedYPosition(x(i) - observer_position(0), 
                                              y(i) - observer_position(1), accu_hazard(i), u(i), 
                                              parameter); 
            recorded_t = t * dt + (y(i) - observer_position(1) - recorded_y) / 
              observer_speed;
            if ((recorded_t < survey_time) & (fabs(x(i) - observer_position(0)) < half_width)) include = true; 
          } else {
            recorded_y = y(i) - observer_position(1); 
            recorded_t = t * dt + GetRecordedT(x(i) - observer_position(0), 
                                               y(i) - observer_position(1), 
                                               accu_hazard(i), 
                                               u(i), 
                                               parameter);
            recorded_r = pow(x(i) - observer_position(0), 2) + recorded_y * recorded_y; 
            if (recorded_r <= half_width * half_width) include = true; 
          }
          if (include) {
            data << transect + 1 << "," << x(i) - observer_position(0) << "," << recorded_y << "," <<
              recorded_t << "\n"; 
          }
        }
        accu_hazard(i) += hazard; 
        // assume uniform initial behaviour distribution 
        if (move == 1) {
          x(i) += sd * sqrt(dt) * R::rnorm(0, 1);
          y(i) += sd * sqrt(dt) * R::rnorm(0, 1);
        }  
      }
      observer_position(1) += observer_speed * dt; 
      arma::uvec outside = find(x < 0 || x > area(0) || y < 0 || y > area(1)); 
      arma::vec u2(arma::randu(N)); 
      u(outside) = u2(outside);
      accu_hazard(outside).fill(0.0);  
      x(find(x < 0)) += area(0); 
      x(find(x > area(0))) -= area(0); 
      y(find(y < 0)) += area(1); 
      y(find(y > area(1))) -= area(1);  
    }  
  }
  data.close(); 
  return 0; 
}