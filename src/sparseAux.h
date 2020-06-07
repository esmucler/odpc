#ifndef SPARSEAUX_H
#define SPARSEAUX_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"

double Positive_Part(const double & x);
double Sign(const double & x);
double Absolute_Value(const double & x);
double Soft_Thresholding(const double & x, const double & gamma);
void Vector_Soft_Thresholding(const double & gamma, arma::vec & a);
#endif