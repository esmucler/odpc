#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"
#include "miscAux.h"

// [[Rcpp::export]]
double getMSE(const arma::mat & resp,
              const arma::mat & Fitted){
  // Get MSE of the reconstruction of resp by Fitted
  // INPUT
  // resp: matrix ot be reconstructed
  // Fitted: matrix of fitted values
  // OUTPUT
  // mse: mean squared error
  int N = resp.n_rows;
  int m = resp.n_cols;
  double mse = accu(pow(resp - Fitted, 2));
  mse /= (N * m);
  return(mse);
}

double getObj(const arma::mat & resp,
              const arma::mat & Fitted,
              const arma::vec & a,
              const double & lambda,
              const double & alpha_en){
  
  int N = resp.n_rows;
  int m = resp.n_cols;
  double obj = getMSE(resp, Fitted) + lambda * alpha_en * norm(a, 1) + lambda * (1-alpha_en) * 0.5 * pow(norm(a, 2),2);
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
                                      const arma::vec & fout,
                                      const double & lambda,
                                      const double & mse){
  arma::field<arma::mat> pre_ret(12, 1);
  arma::vec k1_vec = zeros(1);
  arma::vec k2_vec = zeros(1);
  arma::vec obj_vec = zeros(1);
  arma::vec mse_vec = zeros(1);
  arma::vec criter_vec = zeros(1);
  arma::vec conv_vec = zeros(1);
  k1_vec(0) = k1;
  k2_vec(0) = k2;
  obj_vec(0) = obj;
  mse_vec(0) = mse;
  criter_vec(0) = criter;
  conv_vec(0) = conv;
  
  pre_ret(0, 0) = alpha;
  pre_ret(1, 0) = B;
  pre_ret(2, 0) = k2_vec;
  pre_ret(3, 0) = mse_vec;
  pre_ret(4, 0) = fout;
  pre_ret(5, 0) = res;
  pre_ret(6, 0) = k1_vec;
  pre_ret(7, 0) = criter_vec;
  pre_ret(8, 0) = conv_vec;
  pre_ret(9, 0) = a;
  pre_ret(10, 0) = lambda;
  pre_ret(11, 0) = obj;
  return(pre_ret);
}