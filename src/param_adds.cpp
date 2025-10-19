// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//


// [[Rcpp::export]]
std::complex<double> Compute_charact_funct(double t, double mu, double alpha,
                                           double beta, double sigma){

  const double pi = std::acos(-1);
  double PHI{0};

  if (alpha != 1)
    PHI = tan(pi*alpha/2);
  else
    PHI = -2/pi*log(std::abs(t));

  double angle{0};
  angle = t*mu + pow(std::abs(sigma*t), alpha) * beta * std::abs(t)/t * PHI;
  // t here must not be zero !

  std::complex<double> out{0, 0};
  double cmlx_scale{0};

  cmlx_scale = exp(-pow(std::abs(sigma*t), alpha));
  out = std::polar(cmlx_scale, angle);

  return(out);
}


// [[Rcpp::export]]
std::complex<double> ComplexCF_cpp(double t, Rcpp::NumericVector theta)
{
  double alpha{0};
  double beta{0};
  double sigma{0};
  double mu{0};

  alpha = theta[0];
  beta = theta[1];
  sigma = theta[2];
  mu = theta[3];

  // CheckParametersRange(c(alpha, beta, gamma, delta))
  std::complex<double> out;
  out = Compute_charact_funct(t, mu, alpha, beta, sigma);
  return(out);
}



// [[Rcpp::export]]
arma::mat stbl_param_covmtrx_cpp(double t_par, Rcpp::NumericVector theta){

  std::complex<double> phi_2t;
  std::complex<double> phi_m2t;

  phi_2t = ComplexCF_cpp(2*t_par, theta);
  phi_m2t = ComplexCF_cpp(-2*t_par, theta);

  std::complex<double> phi_t;
  std::complex<double> phi_mt;

  phi_t = ComplexCF_cpp(t_par, theta);
  phi_mt = ComplexCF_cpp(-1*t_par, theta);

  double el_11{0};
  double el_22{0};
  double el_nd{0};
  const std::complex<double> ifour{0, 4};

  el_11 = real(0.25 * (phi_2t + 2.0 + phi_m2t - std::pow((phi_t), 2) - 2.0*phi_t*phi_mt - std::pow((phi_mt), 2)));
  el_22 = real(0.25 * (std::pow((phi_t), 2) - 2.0*phi_t*phi_mt + std::pow((phi_mt), 2)) - 0.25*(phi_2t + phi_m2t - 2.0));
  el_nd = real(1.0/(ifour) * (phi_2t - std::pow((phi_t), 2) - phi_m2t + std::pow((phi_mt), 2)));

  arma::mat mtrx(2, 2, arma::fill::zeros);
  mtrx(0,0) = el_11;
  mtrx(1,1) = el_22;
  mtrx(0,1) = el_nd;
  mtrx(1,0) = el_nd;
  return(mtrx);
}
