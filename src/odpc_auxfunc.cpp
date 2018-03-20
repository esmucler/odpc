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

arma::mat getMatrixXj(const arma::mat & Z, const int & k1, const int & j){
  // Get submatrix Xj made up of rows j+1 to T-k1+j of Z
  int N = Z.n_rows;
  arma::mat Xj = Z.rows(j, N - k1 + j - 1);
  return(Xj);
}

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

arma::mat getMatrixX(const arma::mat & Z, const int & k1){
  // Get X = (Xk1+1, ..., X1)
  int N = Z.n_rows;
  int m = Z.n_cols;
  arma::mat X = zeros(N - k1, m * (k1 + 1));
  for (int h = 0 ; h <= k1; h++) {
    X.cols(h * m , (h + 1) * m - 1) = getMatrixXj(Z, k1, k1 - h);
  }
  return(X);
}


arma::mat getMatrixC_forecast(const arma::mat & Z, const int & k_tot,
                              const int & k1, const int & k2){
  // Get C = [Zk_tot0; ...; Zk10]
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
  int L = f.n_elem; // = T - k1
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

void getMatrixD(const arma::mat & Z2k, const arma::mat & F, arma::mat & outD){
  // Get matrix D of loadings and intercepts. First row contains the intercepts (alpha)
  // INPUT
  // Z2k: submatrix of data matrix Z formed by keeping only rows (k1+k2)+1, ..., T
  // F: current matrix of principal components
  // OUTPUT
  // outD: matrix of loadings and intercepts
  double condi = cond(F);
  if (condi < 1e10){
    outD = solve(F, Z2k);
  } else {
    outD = pinv(F) * Z2k;
  }
}


// [[Rcpp::export]]
arma::mat getFini_forecast(const arma::mat & Z,
                           const int & k1,
                           const int & k2,
                           const int & k_tot) {
  // Get initial matrix F: built using the ordinary principal component with k_1 lags
  int N = Z.n_rows;
  arma::mat Z_cen = Z;
  arma::vec f_ini = zeros(N, 1);
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::rowvec mean_Z = mean(Z);
  Z_cen.each_row() -= mean_Z; 
  svd_econ(U, s, V, Z_cen, "right");
  f_ini = Z * V.col(0);
  arma::mat Fini = zeros(N - k_tot, k2 + 2);
  for (int h = 0; h <= k2; h++){
    Fini.col(h + 1) = f_ini(span(k_tot - h, N - h - 1));
  }
  //First column is filled with ones
  Fini.col(0).fill(1);
  return(Fini);
}

void getVecAini(const arma::mat Z,
                const arma::vec f_ini,
                const int & k1,
                arma::vec & outA){
  int N = Z.n_rows;
  int m = Z.n_cols;
  arma::mat X = zeros(N - k1, m * (k1 + 1));
  X = getMatrixX(Z, k1);
  double condi = cond(X);
  if (condi < 1e10) {
    outA = solve(X, f_ini);
  } else {
    outA = pinv(X) * f_ini;
  }
  outA /= norm(outA);
}


void getVecA(arma::sp_mat & W,
             arma::mat & WC,
             const arma::mat & B,
             const arma::mat & ident,
             const arma::mat & C,
             const arma::mat & vecZ2k,
             arma::vec & outA){
  // Get vector a that defines the principal component
  // INPUT
  // W: reference to W = (B' \kron ident)
  // B: current matrix of loadings
  // ident: (T-(k1+k2)) x (T-(k1+k2))  identity matrix
  // C: = [Z(k1+k2)0; ...; Zk10]
  // vecZ2k: vectorization of Z2k - alphas
  // OUTPUT
  // outA: vector a that defines the principal component
  W = kron(B.t(), ident);
  WC = W * C;
  double condi = cond(WC);
  if (condi < 1e10) {
    outA = solve(WC, vecZ2k);
  } else {
    outA = pinv(WC) * vecZ2k;
  }
  outA /= norm(outA);
}

void getVecAMatD(const arma::mat & Z2k,
                 const arma::mat & matF,
                 const arma::mat & ident,
                 const arma::mat & C,
                 const arma::vec & one,
                 arma::vec & out_WCres,
                 arma::vec & outa,
                 arma::vec & outalpha,
                 arma::mat & outB,
                 arma::mat & outD,
                 arma::vec & vecZ2k,
                 arma::sp_mat & W){
  int k = outD.n_rows - 2;
  getMatrixD(Z2k, matF, outD);
  outB = outD.rows(1, k + 1);
  outalpha = outD.row(0).t();
  vecZ2k = vectorise(Z2k) - kron(outalpha, ident) * one;
  W = kron(outB.t(), ident);
  for (arma::uword iter_a = 0; iter_a < outa.n_rows; iter_a++) {
    outa[iter_a] = 0;
    out_WCres = W * C.col(iter_a);
    outa[iter_a] = dot(out_WCres, vecZ2k - W * (C * outa) );
    outa[iter_a] /=  pow(norm(out_WCres), 2);
  }
  double norma = norm(outa);
  outa /= norma;
  for (arma::uword i = 1; i < outD.n_rows; i++){
    outD.row(i) *= norma; 
  }
}


// [[Rcpp::export]]
double getMSE(const arma::mat & Z2k,
              const arma::mat & Fitted){
  // Get MSE of the reconstruction of Z2k by Fitted
  // INPUT
  // Z2k: submatrix of data matrix Z formed by keeping only rows ()k1+k2)+1, ..., T
  // Fitted: matrix of fitted values
  // OUTPUT
  // mse: mean squared error
  int N = Z2k.n_rows;
  int m = Z2k.n_cols;
  double mse = accu(pow(Z2k - Fitted, 2));
  mse /= (N * m);
  return(mse);
}


// [[Rcpp::export]]
List odpc_priv(const arma::mat & Z,
                const int & k1,
                const int & k2,
                const arma::vec & f_ini,
                const bool & passf_ini,
                const double & tol,
                const int & niter_max,
                const int & method) {
  // This function computes a single ODPC with a given number of lags.
  // INPUT
  // Z: data matrix each column is a different time series
  // k1: number of lags used to define f
  // k2: number of lags used to reconstruct
  // f_ini: initial estimate of f
  // passf_ini: logical: is f_ini being passed?
  // tol: relative precision, stopping criterion
  // niter_max: maximum number of iterations
  // method: 1 =  ALS, 2 = CD in a, LS in B
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
  int k_tot = k1 + k2;
  arma::mat ident = eye(N - k_tot, N - k_tot);
  arma::vec vecZ2k = zeros((N - k_tot) * m, 1); //will store vectorise(Z2k) - alphas
  arma::vec one = ones(N - k_tot, 1);
  arma::vec alpha = zeros(m, 1);
  arma::mat Z2k = zeros(N - k_tot, m);
  arma::mat res = zeros(N - k_tot, m);
  arma::mat C = zeros((N - k_tot) * (k2 + 1), m * (k1 + 1));
  arma::mat B = zeros(k2 + 1, m);
  arma::mat D = zeros(k2 + 2, m);
  arma::mat matF = zeros(N - k_tot, k2 + 2);
  matF.col(0).fill(1); //matF's first column is filled with 0
  arma::vec fout = zeros(N - k1, 1);
  arma::mat Fitted = zeros(N - k_tot, m);
  arma::sp_mat W = sp_mat(m * (N - k_tot), (N - k_tot) * (k2 + 1));
  arma::vec a = zeros(m * (k1 + 1), 1);
  double mse = 0;
  int niter = 0;
  bool conv = false;
  double criter = tol + 1;
  Z2k = getMatrixZj(Z, k_tot, k_tot);
  C = getMatrixC_forecast(Z, k_tot, k1, k2);
  
  // if using ALS method
  if (method == 1){
    arma::mat WC = zeros(m * (N - k_tot), m * (k2 + 1));
    // if initial f was passed, get initial a associated with it
    // else use ordinary PC
    if (passf_ini) {
      getVecAini(Z, f_ini, k1, a);
      getMatrixF(Z, k1, k2, k_tot, a, matF);
      getMatrixD(Z2k, matF, D);
      B = D.rows(1, k2 + 1);
      alpha = D.row(0).t();
      vecZ2k = vectorise(Z2k) - kron(alpha, ident) * one;
    } else {
      matF = getFini_forecast(Z, k1, k2, k_tot);
      getMatrixD(Z2k, matF, D);
      B = D.rows(1, k2 + 1);
      alpha = D.row(0).t();
      vecZ2k = vectorise(Z2k) - kron(alpha, ident) * one;
      getVecA(W, WC, B, ident, C, vecZ2k, a);
      getMatrixF(Z, k1, k2, k_tot, a, matF);
      getMatrixD(Z2k, matF, D);
      B = D.rows(1, k2 + 1);
      alpha = D.row(0).t();
      vecZ2k = vectorise(Z2k) - kron(alpha, ident) * one;
    }
    Fitted = matF * D;
    mse = getMSE(Z2k, Fitted);
    double mse_ini = mse;
    while (niter < niter_max and criter > tol){
      niter += 1;
      getVecA(W, WC, B, ident, C, vecZ2k, a);
      getMatrixF(Z, k1, k2, k_tot, a, matF);
      getMatrixD(Z2k, matF, D);
      B = D.rows(1, k2 + 1);
      alpha = D.row(0).t();
      vecZ2k = vectorise(Z2k) - kron(alpha, ident) * one;
      Fitted = matF * D;
      mse = getMSE(Z2k, Fitted);
      criter = 1 - mse / mse_ini;
      mse_ini = mse;
      // if (niter % 2 == 0){
      //   Rcpp::checkUserInterrupt();
      // }
    }
    // if using mix method
  } else if (method == 2){
    arma::vec out_WCres = zeros((N - k_tot) * m);
    if (passf_ini) {
      getVecAini(Z, f_ini, k1, a);
      getMatrixF(Z, k1, k2, k_tot, a, matF);
    } else {
      matF = getFini_forecast(Z, k1, k2, k_tot);
    }
    mse = 1e10;
    double mse_ini = mse;
    while (niter < niter_max and criter > tol){
      niter += 1;
      getVecAMatD(Z2k, matF, ident, C, one, out_WCres, a, alpha, B, D, vecZ2k, W);
      getMatrixF(Z, k1, k2, k_tot, a, matF);
      Fitted = matF * D;
      mse = getMSE(Z2k, Fitted);
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
    fout(span(k2, N - k1 - 1)) = matF.col(1);
  } else {
    fout = matF.col(1);
  }
  
  res = Z2k - Fitted;
  
  List ret;
  ret["a"] = a;
  ret["alpha"] = alpha;
  ret["B"] = B;
  ret["k1"] = k1;
  ret["k2"] = k2;
  ret["mse"] = mse;
  ret["conv"] = conv;
  ret["criter"] = criter;
  ret["f"] = fout;
  ret["res"] = res;
  return(ret);
}

// [[Rcpp::export]]
List roll_odpc(const arma::field<arma::mat> & data_field,
               const arma::uword & k,
               const int & h,
               const int & window_size,
               const double & tol,
               const int & niter_max,
               const int & method,
               const arma::uword & ncores) {
  // This function computes a ODPC over a rolling window.
  // INPUT
  // data_field: Field with the data to be used for training (could be residuals from previous fits),
  // k: value of k to use (k1==k2)
  // h: steps ahead
  // window_size: the size of the window
  // tol: relative precision, stopping criterion
  // niter_max: maximum number of iterations
  // method: 1 =  ALS, 2 = CD in a, LS in B
  // OUTPUT
  // a list of the same length as windows_size, each entry being an ODPC
  arma::vec nothing = zeros(2);
  List output(window_size);
  # pragma omp parallel for num_threads(ncores)
  for (int ind=0; ind < window_size; ind++){
    output[ind] = odpc_priv(data_field(ind, 0), k, k, nothing, false, tol, niter_max, method);
  }
  return(output);
}

// [[Rcpp::export]]
List grid_odpc(const arma::field<arma::mat> data_field,
               const arma::vec & k_list,
               const arma::uword & h,
               const arma::uword & window_size,
               const double & tol,
               const int & niter_max,
               const int & method,
               const arma::uword & ncores) {
  // This function computes a ODPC for each k in k_list, over a rolling window of size window_size.
  // INPUT
  // data_field: field with the data to be used for training (could be residuals from previous fits),
  // one entry per window size
  // k_list: list of ks to use (k1==k2)
  // h: number of steps ahead, only used for determining the training set
  // window_size: the size of the rolling window
  // tol: relative precision, stopping criterion
  // niter_max: maximum number of iterations
  // method: 1 =  ALS, 2 = CD in a, LS in B
  // OUTPUT
  // a list of the same length as k_list, each entry being a list of ODPCs computed over the rolling window
  // for the corresponding value of k
  
  arma::uword k_list_len = k_list.n_elem;
  List output(k_list_len);

  for (arma::uword ind =0; ind < k_list_len; ind++){
    output[ind] = roll_odpc(data_field, k_list[ind], h, window_size, tol, niter_max, method, ncores);
  }
  return(output);
}