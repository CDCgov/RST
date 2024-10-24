#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec u_update_Z(
  arma::vec& Z,
  const double& sig2,
  const arma::vec& theta,
  const arma::vec& beta,
  const double& tau2,
  const arma::field<arma::uvec> adjacency,
  const arma::vec& num_adj,
  const arma::field<arma::uvec>& island_region,
  const arma::uvec& island_id
) {
  uword num_region = Z.n_elem;
  for (uword reg = 0; reg < num_region; reg++) {
    double var_Z = 1 / (1 / tau2 + num_adj[reg] / sig2);
    double mean_Z = var_Z * ((theta[reg] - beta[island_id[reg]]) / tau2 + sum(Z.elem(adjacency[reg])) / sig2);
    Z[reg] = R::rnorm(mean_Z, sqrt(var_Z));
  }
  Z -= mean(Z);
  return Z;
}

//[[Rcpp::export]]
double u_update_sig2(
  double& sig2,
  const arma::vec& Z,
  const arma::vec& beta,
  const double& tau2,
  const arma::field<arma::uvec> adjacency,
  const arma::vec& num_adj,
  const arma::field<arma::uvec>& island_region,
  const arma::uvec& num_island_region,
  const double& A,
  const double& m0,
  const double& sig_a,
  const double& sig_b,
  const String& method
) {
  uword num_region = Z.n_elem;
  uword num_island = island_region.n_elem;
  double sum_adj = 0;
  for (uword reg = 0; reg < num_region; reg++) {
    sum_adj += Z[reg] * sum(Z.elem(adjacency[reg]));
  }
  double a_sig = (num_region - num_island) / 2 + sig_a;
  double b_sig = 1 / ((sum(pow(Z, 2) % num_adj) - sum_adj) / 2 + sig_b);
  if (A >= 1e2) {
    sig2 = 1 / R::rgamma(a_sig, b_sig);
  } else if (A < 1e2) {
    double sig_thres = 0;
    if (method == "binom") {
      double pi = sum(beta % num_island_region / num_region);
      pi = exp(pi) / (1 + exp(pi));
      sig_thres = (1 / ((A + pi) * (1 - pi)) - tau2 * (1 + 1 / m0)) * m0;
    } else if (method == "pois") {
      sig_thres = (log(1 / A + 1) - tau2 * (1 + 1 / m0)) * m0;
    }
    sig_thres = (sig_thres < 0) ? 0 : sig_thres;
    double u = R::runif(0, R::pgamma(1 / sig_thres, a_sig, b_sig, true, false));
    sig2 = 1 / R::qgamma(u, a_sig, b_sig, true, false);
  } 
  return sig2;
}

//[[Rcpp::export]]
double u_update_tau2(
  double& tau2,
  const arma::vec& theta,
  const arma::vec& beta,
  const arma::vec& Z,
  const double& sig2,
  const arma::uvec& num_island_region,
  const arma::uvec& island_id,
  const double& A,
  const double& m0,
  const double& tau_a,
  const double& tau_b,
  const String& method
) {
  uword num_region = Z.n_elem;
  double a_tau = num_region / 2 + tau_a;
  double b_tau = 1 / (sum(pow(theta - beta.elem(island_id) - Z, 2)) / 2 + tau_b);
  if (A >= 1e2) {
    tau2 = 1 / R::rgamma(a_tau, b_tau);
  } else if (A < 1e2) {
    double tau_thres = 0;
    if (method == "binom") {
      double pi = sum(beta % num_island_region / num_region);
      pi = exp(pi) / (1 + exp(pi));
      tau_thres = (1 / ((A + pi) * (1 - pi)) - sig2 / m0) / (1 + 1 / m0);
    } else if (method == "pois") {
      tau_thres = (log(1 / A + 1) - sig2 / m0) / (1 + 1 / m0);
    }
    tau_thres = (tau_thres < 0) ? 0 : tau_thres;
    double u = R::runif(0, R::pgamma(1 / tau_thres, a_tau, b_tau, true, false));
    tau2 = 1 / R::qgamma(u, a_tau, b_tau, true, false);
  }
  return tau2;
}

//[[Rcpp::export]]
arma::vec u_update_theta(
  arma::vec& theta,
  arma::vec& t_accept,
  const arma::vec& Y,
  const arma::vec& n,
  const arma::vec& Z,
  const arma::vec& beta,
  const double& tau2,
  const arma::vec& theta_sd,
  const arma::uvec& island_id,
  const String& method
) {
  uword num_region = Z.n_elem;
  for (uword reg = 0; reg < num_region; reg++) {
    double theta_star = R::rnorm(theta(reg), theta_sd(reg));
    double rk1 = Y(reg) * (theta_star - theta(reg));
    double rk2 = 0;
    if (method == "binom") {
      rk2 = n(reg) * (log(1 + exp(theta_star)) - log(1 + exp(theta(reg))));
    } else if (method == "pois") {
      rk2 = n(reg) * (exp(theta_star) - exp(theta(reg)));
    }
    double rk3a = pow(theta_star - beta(island_id[reg]) - Z(reg), 2);
    double rk3b = pow(theta(reg) - beta(island_id[reg]) - Z(reg), 2);
    double rk   = exp(rk1 - rk2 - 1 / (2 * tau2) * (rk3a - rk3b));
    if (rk >= R::runif(0, 1)) {
      t_accept(reg)++;
      theta(reg) = theta_star;
    }
  }
  return theta;
}

//[[Rcpp::export]]
arma::vec u_update_beta(
  arma::vec& beta,
  const arma::vec& theta,
  const arma::vec& Z,
  const double& tau2,
  const double& sig2,
  const double& A,
  const double& m0,
  const arma::field<arma::uvec>& island_region,
  const String& method
) {
  uword num_island = island_region.n_elem;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region[isl].n_elem;
    double sd_beta = sqrt(tau2 / num_island_region);
    double mean_beta = mean(theta.elem(island_region[isl]) - Z.elem(island_region[isl]));
    if (A >= 1e2) {
      beta[isl] = R::rnorm(mean_beta, sd_beta);
    } else if (A < 1e2) {
      if (method == "binom") {
        double var_t = tau2 + (tau2 + sig2) / m0;
        double pi_beta = pow(A - 1, 2) + 4 * (A - 1 / var_t);
        double beta_thres = ((1 - A) + sqrt(pi_beta)) / 2;
        beta_thres = log(beta_thres / (1 - beta_thres));
        beta_thres = (beta_thres < 0) ? 0 : beta_thres;
        double beta_max = R::pnorm(beta_thres, mean_beta, sd_beta, true, false);
        if (beta_max > 0) {
          double u = R::runif(0, beta_max);
          beta[isl] = R::qnorm(u, mean_beta, sd_beta, true, false);
        }
      } else if (method == "pois") {
        beta[isl] = R::rnorm(mean_beta, sd_beta);
      }
    }
  }
  return beta;
}