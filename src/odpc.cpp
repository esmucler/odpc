#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"
#include "getMatrices.h"
#include "updateFunctions.h"

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


// [[Rcpp::export]]
arma::field<arma::mat> odpc_priv(const arma::mat & Z,
                                 const arma::mat & resp,
                                 const int & k_tot_max,
                                 const int & k1,
                                 const int & k2,
                                 const arma::uword & num_comp,
                                 const arma::vec & f_ini,
                                 const bool & passf_ini,
                                 const double & tol,
                                 const int & niter_max,
                                 const int & method) {
  // This function computes a single ODPC with a given number of lags.
  // INPUT
  // Z: data matrix each column is a different time series
  // resp: series to be reconstructed; if q components have been computed, this will have
  // N-k_tot_max, where k_tot_max=max(k^i1+k^i2)
  // k1: number of lags used to define f
  // k2: number of lags used to reconstruct
  // num_comp: what component is this?
  // k_tot_max: max(k^i1+k^i2)
  // f_ini: initial estimate of f
  // passf_ini: logical: is f_ini being passed?
  // tol: relative precision, stopping criterion
  // niter_max: maximum number of iterations
  // method: 1 =  ALS, 2 = CD in a, LS in B, 3 = GD in A, LS in B
  // OUTPUT
  // k1: number of lags used to define f
  // k2: number of lags used to reconstruct
  // a: vector to construct the principal component
  // alpha: vector of intercepts corresponding to the principal component
  // B: matrix of loadings corresponding to the principal component
  // mse:  mean squared error
  // conv: logical. Did the iterations converge?
  // res: matrix of residuals
  // f: matrix F
  // criter: last value of 1-mse1/mse0
  int N = Z.n_rows;
  int m = Z.n_cols;
  
  arma::mat ident = eye(N - k_tot_max, N - k_tot_max);
  arma::vec vecresp = zeros((N - k_tot_max) * m, 1); //will store vectorise(resp) - alphas
  arma::vec one = ones(N - k_tot_max, 1);
  arma::vec alpha = zeros(m, 1);
  arma::mat res = zeros(N - k_tot_max, m);
  arma::mat C = zeros((N - k_tot_max) * (k2 + 1), m * (k1 + 1));
  arma::mat B = zeros(k2 + 1, m);
  arma::mat D = zeros(k2 + 2, m);
  arma::mat matF = zeros(N - k_tot_max, k2 + 2);
  matF.col(0).fill(1); //matF's first column is filled with 0
  arma::vec fout = zeros(N - (k_tot_max - k2), 1);
  arma::mat Fitted = zeros(N - k_tot_max, m);
  arma::sp_mat W = sp_mat(m * (N - k_tot_max), (N - k_tot_max) * (k2 + 1));
  arma::vec a = zeros(m * (k1 + 1), 1);
  double mse = 0;
  int niter = 0;
  bool conv = false;
  double criter = tol + 1;
  C = getMatrixC_forecast(Z, k_tot_max, k1, k2);
  // if using ALS method
  if (method == 1){
    arma::mat WC = zeros(m * (N - k_tot_max), m * (k2 + 1));
    matF = getFini_forecast(Z, resp, k1, k2, num_comp);
    getMatrixD(resp, matF, D);
    B = D.rows(1, k2 + 1);
    alpha = D.row(0).t();
    vecresp = vectorise(resp) - kron(alpha, ident) * one;
    getVecA(W, WC, B, ident, C, vecresp, a);
    getMatrixF(Z, k1, k2, k_tot_max, a, matF);
    getMatrixD(resp, matF, D);
    B = D.rows(1, k2 + 1);
    alpha = D.row(0).t();
    vecresp = vectorise(resp) - kron(alpha, ident) * one;
    Fitted = matF * D;
    mse = getMSE(resp, Fitted);
    double mse_ini = mse;
    while (niter < niter_max and criter > tol){
      niter += 1;
      getVecA(W, WC, B, ident, C, vecresp, a);
      getMatrixF(Z, k1, k2, k_tot_max, a, matF);
      getMatrixD(resp, matF, D);
      B = D.rows(1, k2 + 1);
      alpha = D.row(0).t();
      vecresp = vectorise(resp) - kron(alpha, ident) * one;
      Fitted = matF * D;
      mse = getMSE(resp, Fitted);
      criter = 1 - mse / mse_ini;
      mse_ini = mse;
      // if (niter % 2 == 0){
      //   Rcpp::checkUserInterrupt();
      // }
    }
    // if using mix method
  } else if (method == 2){
    arma::vec out_WCres = zeros((N - k_tot_max) * m);
    matF = getFini_forecast(Z, resp, k1, k2, num_comp);
    mse = 1e10;
    double mse_ini = mse;
    while (niter < niter_max and criter > tol){
      niter += 1;
      getVecAMatD(resp, matF, ident, C, one, out_WCres, a, alpha, B, D, vecresp, W);
      getMatrixF(Z, k1, k2, k_tot_max, a, matF);
      Fitted = matF * D;
      mse = getMSE(resp, Fitted);
      criter = 1 - mse / mse_ini;
      mse_ini = mse;
      // if (niter % 2 == 0){
      //   Rcpp::checkUserInterrupt();
      // }
    }
    // if using gradient method
  } else if (method == 3){
    arma::mat WC = zeros(m * (N - k_tot_max), m * (k2 + 1));
    matF = getFini_forecast(Z, resp, k1, k2, num_comp);
    mse = 1e10;
    double mse_ini = mse;
    while (niter < niter_max and criter > tol){
      niter += 1;
      getVecAMatD_grad(resp, matF, ident, C, one, WC, a, alpha, B, D, vecresp, W);
      getMatrixF(Z, k1, k2, k_tot_max, a, matF);
      Fitted = matF * D;
      mse = getMSE(resp, Fitted);
      // cout << mse << '\n';
      criter = 1 - mse / mse_ini;
      mse_ini = mse;
      // if (niter % 2 == 0){
      //   Rcpp::checkUserInterrupt();
      // }
    }
  } 

  // check convergence
  if (niter < niter_max) {
    conv = true;
  }

  if (k2 > 0) {
    fout(span(0, k2 - 1)) = matF(span(0, k2 - 1), k2 + 1);
    fout(span(k2, N - (k_tot_max - k2) - 1)) = matF.col(1);
  } else {
    fout = matF.col(1);
  }
  
  res = resp - Fitted;
  
  arma::vec k1_vec = zeros(1);
  arma::vec k2_vec = zeros(1);
  arma::vec mse_vec = zeros(1);
  arma::vec criter_vec = zeros(1);
  arma::vec conv_vec = zeros(1);
  k1_vec(0) = k1;
  k2_vec(0) = k2;
  mse_vec(0) = mse;
  criter_vec(0) = criter;
  conv_vec(0) = conv;
  
  arma::field<arma::mat> pre_ret(10, 1);
  
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
  return(pre_ret);
}

// [[Rcpp::export]]
arma::field<arma::field<arma::mat> > roll_odpc(const arma::field<arma::mat> & data_field,
                                               const arma::field<arma::mat> & response_field,
                                               const arma::uword & k,
                                               const arma::uword & k_tot_max,
                                               const arma::uword & num_comp,
                                               const int & window_size,
                                               const double & tol,
                                               const int & niter_max,
                                               const int & method,
                                               const arma::uword & ncores) {
  // This function computes a ODPC over a rolling window.
  // INPUT
  // data_field: Field with the data to be used for training
  // response_field: Field with the data to be reconstructed
  // k: value of k to use (k1==k2)
  // k_tot_max: current max value of k1+k2
  // num_comp : what component is this?
  // window_size: the size of the window
  // tol: relative precision, stopping criterion
  // niter_max: maximum number of iterations
  // method: 1 =  ALS, 2 = CD in a, LS in B, 3 = gradient descent in a, LS in B
  // OUTPUT
  // a list of the same length as windows_size, each entry being an ODPC
  arma::vec nothing = zeros(2);
  arma::field<arma::field<arma::mat> >  output(window_size, 1);
  # pragma omp parallel for num_threads(ncores)
  for (int ind=0; ind < window_size; ind++){
    output(ind, 0) = odpc_priv(data_field(ind, 0), response_field(ind, 0), k_tot_max, k, k,
                               num_comp, nothing, false, tol, niter_max, method);
  }
  # pragma omp barrier
  return(output);
}