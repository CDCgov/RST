#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::cube get_regs(const arma::cube& arr, const arma::uvec& ind) {
  arma::cube arr_sub(ind.n_elem, arr.n_cols, arr.n_slices);
  for (uword reg = 0; reg < ind.n_elem; reg++) {
    arr_sub.row(reg) = arr.row(ind[reg]);
  }
  return arr_sub;
}

arma::vec get_grp(const arma::cube& arr, const arma::uword& reg, const arma::uword time) {
  arma::vec arr_sub(arr.n_cols);
  for (uword grp = 0; grp < arr.n_cols; grp++) {
    arr_sub(grp) = arr(reg, grp, time);
  }
  return arr_sub;
}

arma::field<arma::mat> Sig_eta_i(const arma::cube& G, const arma::vec& rho) {
  uword num_group = rho.n_elem;
  uword num_time  = G.n_slices;
  mat r  = repmat(rho, 1, num_group);
  mat sr = sqrt(1 - pow(r, 2));
  field<mat> Sei(num_time, num_time);
  Sei(0, 0) = inv(G.slice(0));
  for (uword time = 1; time < num_time; time++) {
    mat Gi = inv(G.slice(time));
    Sei(time - 1, time - 1) += ( r / sr).t() % (r / sr % Gi);
    Sei(time    , time    )  = ( 1 / sr).t() % (1 / sr % Gi);
    Sei(time - 1, time    )  = (-r / sr)     % (1 / sr % Gi).t();
    Sei(time    , time - 1)  = (-r / sr).t() % (1 / sr % Gi);
  }
  return Sei;
}

arma::field<arma::mat> Sig_eta(const arma::field<arma::mat> Sein) {
  uword num_time = Sein.n_rows;
  field<mat> SeSein(num_time, num_time);
  for (uword time = 0; time < num_time; time++) {
    if (time > 0) {
      SeSein(time, time - 1) = inv(Sein(time, time)) * Sein(time, time - 1);
    }
    if (time < num_time - 1) {
      SeSein(time, time + 1) = inv(Sein(time, time)) * Sein(time, time + 1);
    }
  }
  return SeSein;
}

arma::mat cpp_rmvnorm(const arma::vec& mean, const arma::mat& covar) {
  vec out = mean;
  vec rand = Rcpp::rnorm(covar.n_cols, 0, 1);
  out += covar * rand;
  return out;
}

arma::mat geteig(const arma::mat& covar) {
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, covar);
  eigvec *= eigvec.t() % repmat(sqrt(eigval), 1, covar.n_cols);
  return eigvec.t();
}