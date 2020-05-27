#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"

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

void getMatrixF(const arma::mat & Z, const int & k1,
                const int & k2, const int & k_tot,
                const arma::vec & a, arma::mat & outF){
  // Get matrix F whose columns are 1 and f_{j} for j = k_tot, ..., k1
  for (int h = 0; h <= k2; h++){
    //first column is already filled with ones
    outF.col(h + 1) = getMatrixZj0(Z, k1, k_tot, k_tot - h) * a;
  }
}

void getMatrixD(const arma::mat & resp, const arma::mat & F, arma::mat & outD){
  // Get matrix D of loadings and intercepts. First row contains the intercepts (alpha)
  // INPUT
  // resp: data to be reconstructed
  // F: current matrix of principal components
  // OUTPUT
  // outD: matrix of loadings and intercepts
  double condi = cond(F);
  if (condi < 1e10){
    outD = solve(F, resp);
  } else {
    outD = pinv(F) * resp;
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

void getVecA(arma::sp_mat & W,
             arma::mat & WC,
             const arma::mat & B,
             const arma::mat & ident,
             const arma::mat & C,
             const arma::mat & vecresp,
             arma::vec & outA){
  // Get vector a that defines the principal component
  // INPUT
  // W: reference to W = (B' \kron ident)
  // B: current matrix of loadings
  // ident: (T-(k1+k2)) x (T-(k1+k2))  identity matrix
  // C: = [Z(k1+k2)0; ...; Zk10]
  // vecresp: vectorization of resp - alphas
  // OUTPUT
  // outA: vector a that defines the principal component
  W = kron(B.t(), ident);
  WC = W * C;
  double condi = cond(WC);
  if (condi < 1e10) {
    outA = solve(WC, vecresp);
  } else {
    outA = pinv(WC) * vecresp;
  }
  outA /= norm(outA);
}

void getVecAMatD(const arma::mat & resp,
                 const arma::mat & matF,
                 const arma::mat & ident,
                 const arma::mat & C,
                 const arma::vec & one,
                 arma::vec & out_WCres,
                 arma::vec & outa,
                 arma::vec & outalpha,
                 arma::mat & outB,
                 arma::mat & outD,
                 arma::vec & vecresp,
                 arma::sp_mat & W){
  int k = outD.n_rows - 2;
  getMatrixD(resp, matF, outD);
  outB = outD.rows(1, k + 1);
  outalpha = outD.row(0).t();
  vecresp = vectorise(resp) - kron(outalpha, ident) * one;
  W = kron(outB.t(), ident);
  for (arma::uword iter_a = 0; iter_a < outa.n_rows; iter_a++) {
    outa[iter_a] = 0;
    out_WCres = W * C.col(iter_a);
    outa[iter_a] = dot(out_WCres, vecresp - W * (C * outa) );
    outa[iter_a] /=  pow(norm(out_WCres), 2);
  }
  double norma = norm(outa);
  outa /= norma;
  for (arma::uword i = 1; i < outD.n_rows; i++){
    outD.row(i) *= norma; 
  }
}

void getVecAMatD_grad(const arma::mat & resp,
                     const arma::mat & matF,
                     const arma::mat & ident,
                     const arma::mat & C,
                     const arma::vec & one,
                     arma::mat & out_WC,
                     arma::vec & outa,
                     arma::vec & outalpha,
                     arma::mat & outB,
                     arma::mat & outD,
                     arma::vec & vecresp,
                     arma::sp_mat & W){
  int k = outD.n_rows - 2;
  getMatrixD(resp, matF, outD);
  // outD = outD + 2 * eta * matF.t() * (resp - matF * outD);
  outB = outD.rows(1, k + 1);
  outalpha = outD.row(0).t();
  vecresp = vectorise(resp) - kron(outalpha, ident) * one;
  W = kron(outB.t(), ident);
  out_WC = W * C;
  arma::vec grad = (-2) * out_WC.t() *  (vecresp - out_WC * outa);
  double step = (0.5) * pow(norm(grad), 2)/pow(norm(out_WC * grad), 2);
  outa = outa - step * grad;
  // outa = outa + 2 * eta * C.t() * W.t() * vecresp - C.t() * W.t() * W * C * outa;
  double norma = norm(outa);
  outa /= norma;
  for (arma::uword i = 1; i < outD.n_rows; i++){
    outD.row(i) *= norma;
  }
}


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