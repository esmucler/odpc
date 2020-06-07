#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"
#include "getMatrices.h"


arma::mat getMatrixZj(const arma::mat & Z, const int & k_tot, const int & j){
  // Get submatrix Zj made up of rows j+1 to N-k_tot+j of Z
  int N = Z.n_rows;
  arma::mat Zj = Z.rows(j, N - k_tot + j - 1);
  return(Zj);
}

// [[Rcpp::export]]
arma::mat getMatrixZj0(const arma::mat & Z, const int & k1,
                       const int & k_tot, const int & j){
  // Get Zj0 = (Zj, ..., Zj-k1)
  int N = Z.n_rows;
  int m = Z.n_cols;
  arma::mat Zj0 = zeros(N - k_tot, m * (k1 + 1));
  for (int h = 0 ; h <= k1; h++) {
    Zj0.cols(h * m , (h + 1) * m - 1) = getMatrixZj(Z, k_tot, j - h);
  }
  return(Zj0);
}

arma::mat getMatrixC_forecast(const arma::mat & Z, const int & k_tot,
                              const int & k1, const int & k2){
  // Get C = [Zk_tot;k1 ...; Zktot-k2;k1]
  int N = Z.n_rows;
  int m = Z.n_cols;
  arma::mat C = zeros((N - k_tot) * (k2 + 1), m * (k1 + 1));
  for (int h = 0; h <= k2; h++) {
    C.rows(h * (N - k_tot), (h + 1) * (N - k_tot) - 1) = getMatrixZj0(Z, k1, k_tot, k_tot - h);
  }
  return(C);
}



// [[Rcpp::export]]
arma::mat getMatrixFore(const arma::vec & f,
                        const int & k2,
                        const int & h){
  // Get matrix F whose columns are 1 and f_{j} for j = k_tot, ..., k1
  // used in forecast.odpcs
  int L = f.n_elem;
  arma::mat outF = zeros(h, k2 + 2);
  for (int i = 0; i <= k2; i++){
    outF.col(i + 1) = f(span(L - h - i, L - 1 - i));
  }
  outF.col(0).fill(1);
  return(outF);
}

// [[Rcpp::export]]
arma::mat getMatrixFitted(const arma::vec & f,
                          const int & k1,
                          const int & k2){
  // Get matrix F whose columns are 1 and f_{j} for j = k_tot, ..., k1
  // used in fitted.odpc
  int L = f.n_elem; // = T - (k_tot_max - k2)
  arma::mat outF = zeros(L - k2, k2 + 2);
  for (int i = 0; i <= k2; i++){
    outF.col(i + 1) = f(span(k2 - i, L - 1 - i));
  }
  outF.col(0).fill(1);
  return(outF);
}

void getMatrixF(const arma::mat & Z, const int & k1,
                const int & k2, const int & k_tot,
                const arma::vec & a, arma::mat & outF){
  // Get matrix F whose columns are 1 and f_{j} for j = k_tot, ..., k1
  for (int h = 0; h <= k2; h++){
    //first column is already filled with ones
    outF.col(h + 1) = getMatrixZj0(Z, k1, k_tot, k_tot - h) * a;
  }
}

// [[Rcpp::export]]
arma::mat getFini_forecast(const arma::mat & Z,
                           const arma::mat & resp,
                           const int & k1,
                           const int & k2,
                           const arma::uword & num_comp) {
  // Get initial matrix F: built using the ordinary principal component with k_1 lags
  int L = resp.n_rows; // N - k_tot_max
  int N = Z.n_rows;
  arma::mat Fini = zeros(L, k2 + 2);
  // if computing the first component, use the first principal component of the full data
  if (num_comp==1){
    int k_tot = k1 + k2;
    arma::mat Z_cen = Z;
    arma::vec f_ini = zeros(N, 1);
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::rowvec mean_Z = mean(Z);
    Z_cen.each_row() -= mean_Z; 
    svd_econ(U, s, V, Z_cen, "right");
    f_ini = Z * V.col(0);
    for (int h = 0; h <= k2; h++){
      Fini.col(h + 1) = f_ini(span(k_tot - h, N - h - 1));
    }
  }
  // else, use the first principal component of current response;
  // it has nrow N-k_max, f_ini has N-k_max-k2. Fill the missing k2 entries with zeros
  if (num_comp > 1){
    arma::vec f_ini = zeros(L + k2, 1);
    arma::mat resp_cen = resp;
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::rowvec mean_resp = mean(resp);
    resp_cen.each_row() -= mean_resp; 
    svd_econ(U, s, V, resp_cen, "right");
    f_ini(span(k2, L + k2 - 1)) = resp * V.col(0);
    for (int h = 0; h <= k2; h++){
      Fini.col(h + 1) = f_ini(span(k2 - h, L + k2 - h - 1));
    }
  }
  //First column is filled with ones
  Fini.col(0).fill(1);
  return(Fini);
}

// [[Rcpp::export]]
arma::mat getMatrixF_sparse_forecast(const arma::mat & Z, const int & k1,
                                     const int & k2, const int & k_tot,
                                     const arma::vec & a){
  // Get matrix F whose columns are 1 and f_{j} for j = k_tot, ..., k1.
  // This is just a version of getMatrixF without reference outputs that
  // is used in cv.sparse_odpc to get a component using an old and new data
  int N = Z.n_rows;
  arma::mat outF = zeros(N - k_tot, k2 + 2);
  for (int h = 0; h <= k2; h++){
    //first column is already filled with ones
    outF.col(h + 1) = getMatrixZj0(Z, k1, k_tot, k_tot - h) * a;
  }
  outF.col(0).fill(1);
  return(outF);
}