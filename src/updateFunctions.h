#ifndef UPDATEFUNCTIONS_H
#define UPDATEFUNCTIONS_H
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"
#include "getMatrices.h"

void getMatrixD(const arma::mat & resp, const arma::mat & F, arma::mat & outD);

void getVecA(arma::sp_mat & W,
             arma::mat & WC,
             const arma::mat & B,
             const arma::mat & ident,
             const arma::mat & C,
             const arma::mat & vecresp,
             arma::vec & outA);

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
                 arma::sp_mat & W);

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
                      arma::sp_mat & W);
#endif