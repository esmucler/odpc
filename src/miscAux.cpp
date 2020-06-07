#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"
#include "miscAux.h"

double getObj(const arma::mat & resp,
              const arma::mat & Fitted,
              const arma::vec & a,
              const double & lambda){
  
  int N = resp.n_rows;
  int m = resp.n_cols;
  double obj = accu(pow(resp - Fitted, 2)) / (N * m) + lambda * norm(a, 1);
  return(obj);
}

arma::field<arma::mat> process_output(const int & k1,
                                      const int & k2,
                                      const double & obj,
                                      const double & criter,
                                      const bool & conv,
                                      const arma::vec & alpha,
                                      const arma::vec & a,
                                      const arma::mat & B,
                                      const arma::mat & res,
                                      const arma::vec & fout){
  arma::field<arma::mat> pre_ret(10, 1);
  arma::vec k1_vec = zeros(1);
  arma::vec k2_vec = zeros(1);
  arma::vec obj_vec = zeros(1);
  arma::vec criter_vec = zeros(1);
  arma::vec conv_vec = zeros(1);
  k1_vec(0) = k1;
  k2_vec(0) = k2;
  obj_vec(0) = obj;
  criter_vec(0) = criter;
  conv_vec(0) = conv;
  
  pre_ret(0, 0) = alpha;
  pre_ret(1, 0) = B;
  pre_ret(2, 0) = k2_vec;
  pre_ret(3, 0) = obj_vec;
  pre_ret(4, 0) = fout;
  pre_ret(5, 0) = res;
  pre_ret(6, 0) = k1_vec;
  pre_ret(7, 0) = criter_vec;
  pre_ret(8, 0) = conv_vec;
  pre_ret(9, 0) = a;
  return(pre_ret);
}