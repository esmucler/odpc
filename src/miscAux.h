#ifndef MISCAUX_H
#define MISCAUX_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"

arma::field<arma::mat> process_output(const int & k1,
                                      const int & k2,
                                      const double & obj,
                                      const double & criter,
                                      const bool & conv,
                                      const arma::vec & alpha,
                                      const arma::vec & a,
                                      const arma::mat & B,
                                      const arma::mat & res,
                                      const arma::vec & fout,
                                      const double & lambda,
                                      const double & mse);
double getObj(const arma::mat & resp,
              const arma::mat & Fitted,
              const arma::vec & a,
              const double & lambda,
              const double & alpha_en);
double getMSE(const arma::mat & resp,
              const arma::mat & Fitted);
#endif