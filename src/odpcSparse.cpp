#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"
#include "getMatrices.h"
#include "sparseAux.h"
#include "updateFunctions.h"
#include "miscAux.h"

// Solve sparse odpc for one fixed lambda
void solve_sparse_odpc(const arma::mat & Z,
                       const arma::mat & resp,
                       const double & lambda,
                       const double & alpha_en,
                       const int & k_tot_max,
                       const int & k1,
                       const int & k2,
                       const arma::mat & ident,
                       const arma::mat & C,
                       const arma::mat & one,
                       const double & tol,
                       const int & niter_max,
                       arma::mat & WC,
                       arma::vec & a,
                       arma::vec & alpha,
                       arma::mat & B,
                       arma::mat & D,
                       arma::mat & matF,
                       arma::vec & fout,
                       arma::vec & vecresp,
                       arma::sp_mat & W,
                       arma::mat & Fitted,
                       double & obj,
                       double & criter, 
                       bool & conv,
                       double & mse){
  
  double obj_ini = 1e10;
  int niter = 0;
  criter = tol + 1;
  
  getMatrixF(Z, k1, k2, k_tot_max, a, matF);
  
  while (niter < niter_max and criter > tol){
    niter += 1;
    getVecAMatD_grad(resp, matF, ident, C, one, lambda, alpha_en, WC, a, alpha, B, D, vecresp, W);
    getMatrixF(Z, k1, k2, k_tot_max, a, matF);
    Fitted = matF * D;
    obj = getObj(resp, Fitted, a, lambda, alpha_en);
    criter = 1 - obj / obj_ini;
    obj_ini = obj;
    // if (niter % 2 == 0){
    //   Rcpp::checkUserInterrupt();
    // }
  }
  mse = getMSE(resp, Fitted);
  // check convergence
  if (niter < niter_max) {
    conv = true;
  }
  if (k2 > 0) {
    int N = Z.n_rows;
    fout(span(0, k2 - 1)) = matF(span(0, k2 - 1), k2 + 1);
    fout(span(k2, N - (k_tot_max - k2) - 1)) = matF.col(1);
  } else {
    fout = matF.col(1);
  }
}

// Solve sparse odpc over a grid of lambda values using warm starts
void solve_sparse_odpc_grid(const arma::mat & Z,
                             const arma::mat & resp,
                             const arma::vec & lambda_grid,
                             const double & alpha_en,
                             const int & k_tot_max,
                             const int & k1,
                             const int & k2,
                             const arma::mat & ident,
                             const arma::mat & C,
                             const arma::mat & one,
                             const double & tol,
                             const int & niter_max,
                             arma::mat & WC,
                             arma::vec & a,
                             arma::vec & alpha,
                             arma::mat & B,
                             arma::mat & D,
                             arma::mat & matF,
                             arma::vec & fout,
                             arma::vec & vecresp,
                             arma::sp_mat & W,
                             arma::mat & Fitted,
                             double & obj,
                             double & criter,
                             bool & conv,
                             double & mse,
                             arma::field<arma::field<arma::mat>> & ret){
  
  int num_lambda =  ret.n_elem;
  arma::mat res = zeros(1, 1);
  for (arma::uword h=0; h < num_lambda; h++){
    solve_sparse_odpc(Z, resp, lambda_grid[h], alpha_en, k_tot_max, k1, k2, ident, C, one, tol, niter_max, WC, a, alpha,
                      B, D, matF, fout, vecresp, W, Fitted, obj, criter, conv, mse);
    ret[h] = process_output(k1, k2, obj, criter, conv, alpha, a, B, res, fout, lambda_grid[h], mse);
    
  }
}

// [[Rcpp::export]]
arma::field<arma::field<arma::mat>> sparse_odpc_priv(const arma::mat & Z,
                                        const arma::mat & resp,
                                        const int & k_tot_max,
                                        const int & k1,
                                        const int & k2,
                                        const double & tol,
                                        const double & eps,
                                        const int & niter_max,
                                        const arma::vec & a_ini,
                                        const arma::mat & D_ini,
                                        const int & num_lambda_in,
                                        const bool & pass_grid,
                                        const arma::vec & lambda_grid_in,
                                        const double & alpha_en) {
  // This function computes a single sparse ODPC with a given number of lags.
  // INPUT
  // Z: data matrix each column is a different time series
  // resp: series to be reconstructed; if q components have been computed, this will have
  // N-k_tot_max, where k_tot_max=max(k^i1+k^i2)
  // k1: number of lags used to define f
  // k2: number of lags used to reconstruct
  // k_tot_max: max(k^i1+k^i2)
  // tol: relative precision, stopping criterion
  // niter_max: maximum number of iterations
  // a_ini: starting values for vector to construct the principal component
  // D_ini: starting values for matrix of loadings corresponding to the principal component
  // num_lambda: number of l1 penalty constants
  // OUTPUT
  // k1: number of lags used to define f
  // k2: number of lags used to reconstruct
  // a: vector to construct the principal component
  // alpha: vector of intercepts corresponding to the principal component
  // B: matrix of loadings corresponding to the principal component
  // obj:  objective function
  // conv: logical. Did the iterations converge?
  // res: matrix of residuals
  // f: matrix F
  // criter: last value of 1-obj1/obj0
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
  arma::mat WC = zeros(m * (N - k_tot_max), m * (k2 + 1));
  
  C = getMatrixC_forecast(Z, k_tot_max, k1, k2);
  a = a_ini;
  D = D_ini;
  alpha = D.row(0).t();
  B = D.rows(1, k2 + 1);
  bool conv = false;
  double obj = 0;
  double mse = 0;
  double criter = tol + 1;
  
  if (!pass_grid){
    double lambda_max = 0;
    vecresp = vectorise(resp) - kron(alpha, ident) * one;
    W = kron(B.t(), ident);
    WC = W * C;
    double WC_norm =  pow(norm(WC), 2);
    lambda_max = 5 * 1/alpha_en * max(abs(vecresp.t() * WC)) / WC_norm;
    arma::vec lambda_grid = exp(linspace(log(eps * lambda_max), log(lambda_max), num_lambda_in));
    
    arma::field<arma::field<arma::mat>> ret(num_lambda_in);
    solve_sparse_odpc_grid(Z, resp, lambda_grid, alpha_en, k_tot_max, k1,
                           k2, ident, C, one, tol, niter_max,
                           WC, a, alpha, B, D, matF, fout,
                           vecresp, W, Fitted, obj, criter,
                           conv, mse, ret);
    
    return(ret);
  } else {
    int n_lambda = lambda_grid_in.n_elem;
    arma::field<arma::field<arma::mat>> ret(n_lambda);
    solve_sparse_odpc_grid(Z, resp, lambda_grid_in, alpha_en, k_tot_max, k1,
                           k2, ident, C, one, tol, niter_max,
                           WC, a, alpha, B, D, matF, fout,
                           vecresp, W, Fitted, obj, criter,
                           conv, mse, ret);
    
    return(ret);
  }
}

