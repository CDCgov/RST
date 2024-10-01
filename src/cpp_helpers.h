// cpp_helpers.h
#ifndef CPP_HELPERS_H
#define CPP_HELPERS_H

#include <RcppArmadillo.h>

// Use Rcpp and Armadillo namespaces
using namespace Rcpp;
using namespace arma;

// Function declarations

//[[Rcpp::export]]
arma::cube get_regs(const arma::cube& arr, const arma::uvec& ind);

//[[Rcpp::export]]
arma::vec get_grp(const arma::cube& arr, const arma::uword& reg, const arma::uword time);

//[[Rcpp::export]]
arma::field<arma::mat> Sig_eta_i(const arma::cube& G, const arma::vec& rho);

//[[Rcpp::export]]
arma::field<arma::mat> Sig_eta(const arma::field<arma::mat> Sein);

//[[Rcpp::export]]
arma::mat cpp_rmvnorm(const arma::vec& mean, const arma::mat& covar);

//[[Rcpp::export]]
arma::mat geteig(const arma::mat& covar);

#endif // CPP_HELPERS_H
