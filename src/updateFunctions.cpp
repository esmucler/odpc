#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"
#include "getMatrices.h"
#include "updateFunctions.h"
#include "sparseAux.h"

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
                      const double & lambda,
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
  double L = vecresp.n_elem;
  arma::vec grad = (-2) * (1/L) * out_WC.t() *  (vecresp - out_WC * outa);
  double step = (0.5 * L) * pow(norm(grad), 2)/pow(norm(out_WC * grad), 2);
  outa = outa - step * grad;
  // outa = outa + 2 * eta * C.t() * W.t() * vecresp - C.t() * W.t() * W * C * outa;
  if (lambda < 0){
    double norma = norm(outa);
    outa /= norma;
    for (arma::uword i = 1; i < outD.n_rows; i++){
      outD.row(i) *= norma;
    }
  } else {
    double alpha = 0.95;
    Vector_Soft_Thresholding(lambda * step * alpha, outa);
    outa = outa * 1/(1+lambda * (1-alpha));
  }
}