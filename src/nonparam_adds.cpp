// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
double T_sq_nonpar_precomp_cpp(arma::vec means_1, arma::vec means_2,
                               double nrow1, double nrow2,
                               arma::mat cov1, arma::mat cov2) {

  arma::vec diff = means_1 - means_2;
  arma::mat S = (1/(nrow1 + nrow2 - 2)) * (((nrow1 - 1) * cov1) + ((nrow2 - 1) * cov2));
  double res = arma::as_scalar(diff.t() * S.i() * diff);
  res = nrow1*nrow2/(nrow1+nrow2) * res;
  return(res);
}
