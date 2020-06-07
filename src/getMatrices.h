#ifndef GETMATRICES_H
#define GETMATRICES_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"

arma::mat getMatrixZj(const arma::mat & Z, const int & k_tot, const int & j);

arma::mat getMatrixZj0(const arma::mat & Z, const int & k1,
                       const int & k_tot, const int & j);
arma::mat getMatrixC_forecast(const arma::mat & Z, const int & k_tot,
                              const int & k1, const int & k2);
arma::mat getMatrixFore(const arma::vec & f,
                        const int & k2,
                        const int & h);
arma::mat getMatrixFitted(const arma::vec & f,
                          const int & k1,
                          const int & k2);
void getMatrixF(const arma::mat & Z, const int & k1,
                const int & k2, const int & k_tot,
                const arma::vec & a, arma::mat & outF);
arma::mat getFini_forecast(const arma::mat & Z,
                           const arma::mat & resp,
                           const int & k1,
                           const int & k2,
                           const arma::uword & num_comp);
arma::mat getMatrixF_sparse_forecast(const arma::mat & Z, const int & k1,
                                     const int & k2, const int & k_tot,
                                     const arma::vec & a);
#endif