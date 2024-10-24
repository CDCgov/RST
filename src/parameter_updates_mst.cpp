#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::cube mst_update_beta(
  arma::cube& beta,
  const arma::cube& theta,
  const arma::cube& Z,
  const arma::vec& tau2,
  const arma::field<arma::uvec>& island_region
) {
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region[isl].n_elem;
    mat var_beta = repmat(sqrt(tau2 / num_island_region), 1, num_time);
    mat mean_beta = mean(get_regs(theta, island_region[isl]) - get_regs(Z, island_region[isl]), 0);
    beta.row(isl) = mat(num_group, num_time, fill::randn) % var_beta + mean_beta;
  }
  return beta;
}

//[[Rcpp::export]]
arma::cube mst_update_Z(
  arma::cube& Z,
  const arma::cube& G,
  const arma::cube& theta,
  const arma::cube& beta,
  const arma::vec& rho,
  const arma::vec& tau2,
  const arma::field<arma::uvec>& adjacency,
  const arma::vec& num_adj,
  const arma::field<arma::uvec>& island_region,
  const arma::uvec& island_id
) {
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  uword num_island = island_region.n_elem;
  field<mat> Sein = Sig_eta_i(G, rho);
  field<mat> SeSein = Sig_eta(Sein);
  field<mat> Z_cov(num_time, max(num_adj) + 1);
  field<mat> Z_coveig(num_time, max(num_adj) + 1);
  vec unique_num_adj = unique(num_adj);
  for (uword time = 0; time < num_time; time++) {
    for (uword count : unique_num_adj) {
      Z_cov   (time, count) = inv(diagmat(1 / tau2) + count * Sein(time, time));
      Z_coveig(time, count) = geteig(Z_cov(time, count));
    }
  }
  cube rate_diff = theta - get_regs(beta, island_id);
  for (uword reg = 0; reg < num_region; reg++) {
    mat nZm = mean(get_regs(Z, adjacency[reg]), 0);
    for (uword time = 0; time < num_time; time++) {
      vec muZp = nZm.col(time);
      if (time > 0) {
        muZp += SeSein(time, time - 1) * (nZm.col(time - 1) - get_grp(Z, reg, time - 1));
      }
      if (time < num_time - 1) {
        muZp += SeSein(time, time + 1) * (nZm.col(time + 1) - get_grp(Z, reg, time + 1));
      }
      mat Z_mean = Z_cov(time, num_adj[reg]) * (get_grp(rate_diff, reg, time) / tau2 + (num_adj[reg] * Sein(time, time) * muZp));
      vec Z_new = cpp_rmvnorm(Z_mean, Z_coveig(time, num_adj[reg]));
      for (uword grp = 0; grp < num_group; grp++) {
        Z(reg, grp, time) = Z_new(grp);
      }
    }
  }
  cube Zkt(num_island, num_group, num_time);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, island_region[isl]), 0);
  }
  Z -= get_regs(Zkt, island_id);
  return Z;
}

//[[Rcpp::export]]
arma::cube mst_update_G(
    arma::cube& G,
    const arma::cube& Z,
    const arma::mat& Ag,
    const arma::vec& rho,
    const double& G_df,
    const arma::field<arma::uvec>& adjacency,
    const arma::uword& num_island
) {
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  cube Ags(num_group, num_group, num_time, fill::zeros);
  Ags.each_slice() += Ag;
  vec r  = rho;
  vec sr = sqrt(1 - pow(rho, 2));
  for (uword reg = 0; reg < num_region; reg++) {
    double num_adj = adjacency[reg].n_elem;
    mat Zmikt = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
    rowvec Zt = get_grp(Z, reg, 0).t();
    Ags.slice(0) += num_adj * Zmikt.col(0) * Zt;
    for (uword time = 1; time < num_time; time++) {
      vec Zt  = 1 / sr % get_grp(Z, reg, time);
      vec Ztl = r / sr % get_grp(Z, reg, time - 1);
      vec Zm  = 1 / sr % Zmikt.col(time);
      vec Zml = r / sr % Zmikt.col(time - 1);
      Ags.slice(time) += num_adj * ((Zm - Zml) * (Zt - Ztl).t());
    }
  }
  for (uword time = 0; time < num_time; time++) {
    G.slice(time) = riwish((num_region - num_island) + G_df, Ags.slice(time));
  }
  return G;
}

//[[Rcpp::export]]
arma::mat mst_update_Ag(
  arma::mat& Ag,
  const arma::cube& G,
  const arma::mat& Ag_scale,
  const double& G_df,
  const double& Ag_df
) {
  uword num_time = G.n_slices;
  uword num_group = G.n_rows;
  mat Ag_covar(num_group, num_group, fill::zeros);
  Ag_covar += inv(Ag_scale);
  for (uword time = 0; time < num_time; time++) {
    Ag_covar += inv(G.slice(time));
  }
  Ag = rwish(num_time * G_df + Ag_df, inv(Ag_covar));
  return Ag;
}

//[[Rcpp::export]]
arma::vec mst_update_tau2(
  arma::vec& tau2,
  const arma::cube& theta,
  const arma::cube& beta,
  const arma::cube& Z,
  const double& tau_a,
  const double& tau_b,
  const arma::uvec& island_id
) {
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  double tau_shape = num_region * num_time / 2 + tau_a;
  cube square_resid = pow(theta - get_regs(beta, island_id) - Z, 2) / 2;
  for (uword grp = 0; grp < num_group; grp++) {
    double tau_scale = 1 / (accu(square_resid.col(grp)) + tau_b);
    tau2[grp] = 1 / R::rgamma(tau_shape, tau_scale);
  }
  return tau2;
}

//[[Rcpp::export]]
arma::cube mst_update_theta(
  arma::cube& theta,
  arma::cube& t_accept,
  const arma::cube& Y,
  const arma::cube& n,
  const arma::cube& Z,
  const arma::cube& beta,
  const arma::vec& tau2,
  const arma::cube& theta_sd,
  const arma::uvec& island_id,
  const String& method
) {
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  for (uword time = 0; time < num_time; time++) {
    for (uword reg = 0; reg < num_region; reg++) {
      for (uword grp = 0; grp < num_group; grp++) {
        double theta_star = R::rnorm(theta(reg, grp, time), theta_sd(reg, grp, time));
        double rk1 = Y(reg, grp, time) * (theta_star - theta(reg, grp, time));
        double rk2 = 0;
        if (method == "binom") {
          rk2 = n(reg, grp, time) * (log(1 + exp(theta_star)) - log(1 + exp(theta(reg, grp, time))));
        } 
        if (method == "pois") {
          rk2 = n(reg, grp, time) * (exp(theta_star) - exp(theta(reg, grp, time)));
        }
        double rk3a = pow(theta_star            - beta(island_id[reg], grp, time) - Z(reg, grp, time), 2);
        double rk3b = pow(theta(reg, grp, time) - beta(island_id[reg], grp, time) - Z(reg, grp, time), 2);
        double rk   = exp(rk1 - rk2 - 1 / (2 * tau2[grp]) * (rk3a - rk3b));
        if (rk >= R::runif(0, 1)) {
          t_accept(reg, grp, time)++;
          theta(reg, grp, time) = theta_star;
        }
      }
    }
  }
  return theta;
}

//[[Rcpp::export]]
arma::vec mst_update_rho(
  arma::vec& rho,
  arma::vec& r_accept,
  const arma::cube& G,
  const arma::cube& Z,
  const double& rho_a,
  const double& rho_b,
  const arma::vec& rho_sd,
  const arma::field<arma::uvec>& adjacency,
  const arma::uword& num_island
) {
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  vec logit_rho = log(rho / (1 - rho));
  vec rand = Rcpp::rnorm(num_group, 0, 1);
  vec expit_rho = rand % rho_sd + logit_rho;
  vec rho_star_0 = 1 / (1 + exp(-expit_rho));
  vec r(num_group, fill::zeros);
  vec ra(num_group, fill::zeros);
  vec rb(num_group, fill::zeros);
  vec rc(num_group, fill::zeros);
  cube Zm(num_region, num_group, num_time);
  for (uword reg = 0; reg < num_region; reg++) {
    Zm.row(reg) = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
  }
  for (uword grp = 0; grp < num_group; grp++) {
    vec rho_star = rho;
    rho_star[grp] = rho_star_0[grp];
    ra[grp] = (1 - pow(rho[grp], 2)) / (1 - pow(rho_star[grp], 2));
    field<mat> Sein_rho = Sig_eta_i(G, rho);
    field<mat> Sein_rho_star = Sig_eta_i(G, rho_star);
    field<mat> Sein_diff(num_time, num_time);    
    for (uword time1 = 0; time1 < num_time; time1++) {
      uword time1_l = (time1 == 0) ? 0 : time1 - 1;
      uword time2_u = (time1 == num_time - 1) ? num_time : time1 + 2;
      for (uword time2 = time1_l; time2 < time2_u; time2++) {
        Sein_diff(time2, time1) = Sein_rho_star(time2, time1) - Sein_rho(time2, time1);
      }
    }
    for (uword reg = 0; reg < num_region; reg++) {
      for (uword time1 = 0; time1 < num_time; time1++) {
        uword time2_l = (time1 == 0) ? 0 : time1 - 1;
        uword time2_u = (time1 == num_time - 1) ? num_time : time1 + 2;
        for (uword time2 = time2_l; time2 < time2_u; time2++) {
          mat rb_ik = get_grp(Z, reg, time2).t() * Sein_diff(time2, time1) * get_grp(Zm, reg, time1);
          rb[grp] += adjacency[reg].n_elem * rb_ik[0] / 2;
        }
      }
    }
    rc[grp] = pow(rho_star[grp] / rho[grp], rho_a) * pow((1 - rho_star[grp]) / (1 - rho[grp]), rho_b);
    r[grp] = exp((num_region - num_island) * (num_time - 1) / 2 * log(ra[grp]) - rb[grp]) * rc[grp];
    if (r[grp] >= R::runif(0, 1)) {
      r_accept[grp]++;
      rho[grp] = rho_star[grp];
    }
  }
  return rho;
}