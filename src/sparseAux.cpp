#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "sparseAux.h"

// Function to return the positive part of any number
double Positive_Part(const double & x){
  
  double positive = x;
  if (x <= 0){
    positive = 0; 
  }
  return positive;
}

// Function to return the sign of any number - returns a sign of 0 if the numerical argument is 0
double Sign(const double & x){
  
  double sign = 0;
  if (x < 0){
    sign = -1; 
  } 
  if (x > 0){
    sign = 1;
  }
  return sign;
}

// Function that returns the absolute value of any number
double Absolute_Value(const double & x){
  
  double abs_value = x * Sign(x);
  return abs_value;
}

// Function that returns the numerical value from the soft-thresholding operator (takes 2 numerical values as input)
double Soft_Thresholding(const double & x, 
                         const double & gamma){
  
  double soft = 0;
  soft = Sign(x) * Positive_Part(Absolute_Value(x) - gamma);
  return soft;
}

void Vector_Soft_Thresholding(const double & gamma,
                              arma::vec & a){
  
  int n = a.n_elem;
  for (arma::uword h = 0; h < n; h++){
    a[h] = Soft_Thresholding(a[h], gamma);
  }
}