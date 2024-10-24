// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_regs
arma::cube get_regs(const arma::cube& arr, const arma::uvec& ind);
RcppExport SEXP _RSTr_get_regs(SEXP arrSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type arr(arrSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(get_regs(arr, ind));
    return rcpp_result_gen;
END_RCPP
}
// get_grp
arma::vec get_grp(const arma::cube& arr, const arma::uword& reg, const arma::uword time);
RcppExport SEXP _RSTr_get_grp(SEXP arrSEXP, SEXP regSEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type arr(arrSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type reg(regSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_grp(arr, reg, time));
    return rcpp_result_gen;
END_RCPP
}
// Sig_eta_i
arma::field<arma::mat> Sig_eta_i(const arma::cube& G, const arma::vec& rho);
RcppExport SEXP _RSTr_Sig_eta_i(SEXP GSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(Sig_eta_i(G, rho));
    return rcpp_result_gen;
END_RCPP
}
// Sig_eta
arma::field<arma::mat> Sig_eta(const arma::field<arma::mat> Sein);
RcppExport SEXP _RSTr_Sig_eta(SEXP SeinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::mat> >::type Sein(SeinSEXP);
    rcpp_result_gen = Rcpp::wrap(Sig_eta(Sein));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rmvnorm
arma::mat cpp_rmvnorm(const arma::vec& mean, const arma::mat& covar);
RcppExport SEXP _RSTr_cpp_rmvnorm(SEXP meanSEXP, SEXP covarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covar(covarSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rmvnorm(mean, covar));
    return rcpp_result_gen;
END_RCPP
}
// geteig
arma::mat geteig(const arma::mat& covar);
RcppExport SEXP _RSTr_geteig(SEXP covarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type covar(covarSEXP);
    rcpp_result_gen = Rcpp::wrap(geteig(covar));
    return rcpp_result_gen;
END_RCPP
}
// m_update_beta
arma::mat m_update_beta(arma::mat& beta, const arma::mat& theta, const arma::mat& Z, const arma::rowvec& tau2, const arma::field<arma::uvec>& island_region);
RcppExport SEXP _RSTr_m_update_beta(SEXP betaSEXP, SEXP thetaSEXP, SEXP ZSEXP, SEXP tau2SEXP, SEXP island_regionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type island_region(island_regionSEXP);
    rcpp_result_gen = Rcpp::wrap(m_update_beta(beta, theta, Z, tau2, island_region));
    return rcpp_result_gen;
END_RCPP
}
// m_update_Z
arma::mat m_update_Z(arma::mat& Z, const arma::mat& G, const arma::mat& theta, const arma::mat& beta, const arma::vec& tau2, const arma::field<arma::uvec> adjacency, const arma::vec& num_adj, const arma::field<arma::uvec> island_region, const arma::uvec& island_id);
RcppExport SEXP _RSTr_m_update_Z(SEXP ZSEXP, SEXP GSEXP, SEXP thetaSEXP, SEXP betaSEXP, SEXP tau2SEXP, SEXP adjacencySEXP, SEXP num_adjSEXP, SEXP island_regionSEXP, SEXP island_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec> >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type num_adj(num_adjSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec> >::type island_region(island_regionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    rcpp_result_gen = Rcpp::wrap(m_update_Z(Z, G, theta, beta, tau2, adjacency, num_adj, island_region, island_id));
    return rcpp_result_gen;
END_RCPP
}
// m_update_G
arma::mat m_update_G(arma::mat& G, const arma::mat& Z, const double& G_df, const arma::mat& G_scale, const arma::field<arma::uvec>& adjacency, const arma::uword& num_island);
RcppExport SEXP _RSTr_m_update_G(SEXP GSEXP, SEXP ZSEXP, SEXP G_dfSEXP, SEXP G_scaleSEXP, SEXP adjacencySEXP, SEXP num_islandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type G_df(G_dfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type G_scale(G_scaleSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type num_island(num_islandSEXP);
    rcpp_result_gen = Rcpp::wrap(m_update_G(G, Z, G_df, G_scale, adjacency, num_island));
    return rcpp_result_gen;
END_RCPP
}
// m_update_tau2
arma::vec m_update_tau2(arma::vec& tau2, const arma::mat& theta, const arma::mat& beta, const arma::mat& Z, const double& tau_a, const double& tau_b, const arma::uvec& island_id);
RcppExport SEXP _RSTr_m_update_tau2(SEXP tau2SEXP, SEXP thetaSEXP, SEXP betaSEXP, SEXP ZSEXP, SEXP tau_aSEXP, SEXP tau_bSEXP, SEXP island_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_a(tau_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_b(tau_bSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    rcpp_result_gen = Rcpp::wrap(m_update_tau2(tau2, theta, beta, Z, tau_a, tau_b, island_id));
    return rcpp_result_gen;
END_RCPP
}
// m_update_theta
arma::mat m_update_theta(arma::mat& theta, arma::mat& t_accept, const arma::mat& Y, const arma::mat& n, const arma::mat& Z, const arma::mat& beta, const arma::vec& tau2, const arma::mat& theta_sd, const arma::uvec& island_id, const String& method);
RcppExport SEXP _RSTr_m_update_theta(SEXP thetaSEXP, SEXP t_acceptSEXP, SEXP YSEXP, SEXP nSEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP tau2SEXP, SEXP theta_sdSEXP, SEXP island_idSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t_accept(t_acceptSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta_sd(theta_sdSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    Rcpp::traits::input_parameter< const String& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(m_update_theta(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method));
    return rcpp_result_gen;
END_RCPP
}
// mst_update_beta
arma::cube mst_update_beta(arma::cube& beta, const arma::cube& theta, const arma::cube& Z, const arma::vec& tau2, const arma::field<arma::uvec>& island_region);
RcppExport SEXP _RSTr_mst_update_beta(SEXP betaSEXP, SEXP thetaSEXP, SEXP ZSEXP, SEXP tau2SEXP, SEXP island_regionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type island_region(island_regionSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_update_beta(beta, theta, Z, tau2, island_region));
    return rcpp_result_gen;
END_RCPP
}
// mst_update_Z
arma::cube mst_update_Z(arma::cube& Z, const arma::cube& G, const arma::cube& theta, const arma::cube& beta, const arma::vec& rho, const arma::vec& tau2, const arma::field<arma::uvec>& adjacency, const arma::vec& num_adj, const arma::field<arma::uvec>& island_region, const arma::uvec& island_id);
RcppExport SEXP _RSTr_mst_update_Z(SEXP ZSEXP, SEXP GSEXP, SEXP thetaSEXP, SEXP betaSEXP, SEXP rhoSEXP, SEXP tau2SEXP, SEXP adjacencySEXP, SEXP num_adjSEXP, SEXP island_regionSEXP, SEXP island_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type num_adj(num_adjSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type island_region(island_regionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_update_Z(Z, G, theta, beta, rho, tau2, adjacency, num_adj, island_region, island_id));
    return rcpp_result_gen;
END_RCPP
}
// mst_update_G
arma::cube mst_update_G(arma::cube& G, const arma::cube& Z, const arma::mat& Ag, const arma::vec& rho, const double& G_df, const arma::field<arma::uvec>& adjacency, const arma::uword& num_island);
RcppExport SEXP _RSTr_mst_update_G(SEXP GSEXP, SEXP ZSEXP, SEXP AgSEXP, SEXP rhoSEXP, SEXP G_dfSEXP, SEXP adjacencySEXP, SEXP num_islandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Ag(AgSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double& >::type G_df(G_dfSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type num_island(num_islandSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_update_G(G, Z, Ag, rho, G_df, adjacency, num_island));
    return rcpp_result_gen;
END_RCPP
}
// mst_update_Ag
arma::mat mst_update_Ag(arma::mat& Ag, const arma::cube& G, const arma::mat& Ag_scale, const double& G_df, const double& Ag_df);
RcppExport SEXP _RSTr_mst_update_Ag(SEXP AgSEXP, SEXP GSEXP, SEXP Ag_scaleSEXP, SEXP G_dfSEXP, SEXP Ag_dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Ag(AgSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Ag_scale(Ag_scaleSEXP);
    Rcpp::traits::input_parameter< const double& >::type G_df(G_dfSEXP);
    Rcpp::traits::input_parameter< const double& >::type Ag_df(Ag_dfSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_update_Ag(Ag, G, Ag_scale, G_df, Ag_df));
    return rcpp_result_gen;
END_RCPP
}
// mst_update_tau2
arma::vec mst_update_tau2(arma::vec& tau2, const arma::cube& theta, const arma::cube& beta, const arma::cube& Z, const double& tau_a, const double& tau_b, const arma::uvec& island_id);
RcppExport SEXP _RSTr_mst_update_tau2(SEXP tau2SEXP, SEXP thetaSEXP, SEXP betaSEXP, SEXP ZSEXP, SEXP tau_aSEXP, SEXP tau_bSEXP, SEXP island_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_a(tau_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_b(tau_bSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_update_tau2(tau2, theta, beta, Z, tau_a, tau_b, island_id));
    return rcpp_result_gen;
END_RCPP
}
// mst_update_theta
arma::cube mst_update_theta(arma::cube& theta, arma::cube& t_accept, const arma::cube& Y, const arma::cube& n, const arma::cube& Z, const arma::cube& beta, const arma::vec& tau2, const arma::cube& theta_sd, const arma::uvec& island_id, const String& method);
RcppExport SEXP _RSTr_mst_update_theta(SEXP thetaSEXP, SEXP t_acceptSEXP, SEXP YSEXP, SEXP nSEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP tau2SEXP, SEXP theta_sdSEXP, SEXP island_idSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_accept(t_acceptSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type theta_sd(theta_sdSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    Rcpp::traits::input_parameter< const String& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_update_theta(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method));
    return rcpp_result_gen;
END_RCPP
}
// mst_update_rho
arma::vec mst_update_rho(arma::vec& rho, arma::vec& r_accept, const arma::cube& G, const arma::cube& Z, const double& rho_a, const double& rho_b, const arma::vec& rho_sd, const arma::field<arma::uvec>& adjacency, const arma::uword& num_island);
RcppExport SEXP _RSTr_mst_update_rho(SEXP rhoSEXP, SEXP r_acceptSEXP, SEXP GSEXP, SEXP ZSEXP, SEXP rho_aSEXP, SEXP rho_bSEXP, SEXP rho_sdSEXP, SEXP adjacencySEXP, SEXP num_islandSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type r_accept(r_acceptSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho_a(rho_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho_b(rho_bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rho_sd(rho_sdSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type num_island(num_islandSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_update_rho(rho, r_accept, G, Z, rho_a, rho_b, rho_sd, adjacency, num_island));
    return rcpp_result_gen;
END_RCPP
}
// u_update_Z
arma::vec u_update_Z(arma::vec& Z, const double& sig2, const arma::vec& theta, const arma::vec& beta, const double& tau2, const arma::field<arma::uvec> adjacency, const arma::vec& num_adj, const arma::field<arma::uvec>& island_region, const arma::uvec& island_id);
RcppExport SEXP _RSTr_u_update_Z(SEXP ZSEXP, SEXP sig2SEXP, SEXP thetaSEXP, SEXP betaSEXP, SEXP tau2SEXP, SEXP adjacencySEXP, SEXP num_adjSEXP, SEXP island_regionSEXP, SEXP island_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec> >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type num_adj(num_adjSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type island_region(island_regionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    rcpp_result_gen = Rcpp::wrap(u_update_Z(Z, sig2, theta, beta, tau2, adjacency, num_adj, island_region, island_id));
    return rcpp_result_gen;
END_RCPP
}
// u_update_sig2
double u_update_sig2(double& sig2, const arma::vec& Z, const arma::vec& beta, const double& tau2, const arma::field<arma::uvec> adjacency, const arma::vec& num_adj, const arma::field<arma::uvec>& island_region, const arma::uvec& num_island_region, const double& A, const double& m0, const double& sig_a, const double& sig_b, const String& method);
RcppExport SEXP _RSTr_u_update_sig2(SEXP sig2SEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP tau2SEXP, SEXP adjacencySEXP, SEXP num_adjSEXP, SEXP island_regionSEXP, SEXP num_island_regionSEXP, SEXP ASEXP, SEXP m0SEXP, SEXP sig_aSEXP, SEXP sig_bSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec> >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type num_adj(num_adjSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type island_region(island_regionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type num_island_region(num_island_regionSEXP);
    Rcpp::traits::input_parameter< const double& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const double& >::type sig_a(sig_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig_b(sig_bSEXP);
    Rcpp::traits::input_parameter< const String& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(u_update_sig2(sig2, Z, beta, tau2, adjacency, num_adj, island_region, num_island_region, A, m0, sig_a, sig_b, method));
    return rcpp_result_gen;
END_RCPP
}
// u_update_tau2
double u_update_tau2(double& tau2, const arma::vec& theta, const arma::vec& beta, const arma::vec& Z, const double& sig2, const arma::uvec& num_island_region, const arma::uvec& island_id, const double& A, const double& m0, const double& tau_a, const double& tau_b, const String& method);
RcppExport SEXP _RSTr_u_update_tau2(SEXP tau2SEXP, SEXP thetaSEXP, SEXP betaSEXP, SEXP ZSEXP, SEXP sig2SEXP, SEXP num_island_regionSEXP, SEXP island_idSEXP, SEXP ASEXP, SEXP m0SEXP, SEXP tau_aSEXP, SEXP tau_bSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type num_island_region(num_island_regionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    Rcpp::traits::input_parameter< const double& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_a(tau_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_b(tau_bSEXP);
    Rcpp::traits::input_parameter< const String& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(u_update_tau2(tau2, theta, beta, Z, sig2, num_island_region, island_id, A, m0, tau_a, tau_b, method));
    return rcpp_result_gen;
END_RCPP
}
// u_update_theta
arma::vec u_update_theta(arma::vec& theta, arma::vec& t_accept, const arma::vec& Y, const arma::vec& n, const arma::vec& Z, const arma::vec& beta, const double& tau2, const arma::vec& theta_sd, const arma::uvec& island_id, const String& method);
RcppExport SEXP _RSTr_u_update_theta(SEXP thetaSEXP, SEXP t_acceptSEXP, SEXP YSEXP, SEXP nSEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP tau2SEXP, SEXP theta_sdSEXP, SEXP island_idSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type t_accept(t_acceptSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta_sd(theta_sdSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type island_id(island_idSEXP);
    Rcpp::traits::input_parameter< const String& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(u_update_theta(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method));
    return rcpp_result_gen;
END_RCPP
}
// u_update_beta
arma::vec u_update_beta(arma::vec& beta, const arma::vec& theta, const arma::vec& Z, const double& tau2, const double& sig2, const double& A, const double& m0, const arma::field<arma::uvec>& island_region, const String& method);
RcppExport SEXP _RSTr_u_update_beta(SEXP betaSEXP, SEXP thetaSEXP, SEXP ZSEXP, SEXP tau2SEXP, SEXP sig2SEXP, SEXP ASEXP, SEXP m0SEXP, SEXP island_regionSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< const double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const double& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double& >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type island_region(island_regionSEXP);
    Rcpp::traits::input_parameter< const String& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(u_update_beta(beta, theta, Z, tau2, sig2, A, m0, island_region, method));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RSTr_get_regs", (DL_FUNC) &_RSTr_get_regs, 2},
    {"_RSTr_get_grp", (DL_FUNC) &_RSTr_get_grp, 3},
    {"_RSTr_Sig_eta_i", (DL_FUNC) &_RSTr_Sig_eta_i, 2},
    {"_RSTr_Sig_eta", (DL_FUNC) &_RSTr_Sig_eta, 1},
    {"_RSTr_cpp_rmvnorm", (DL_FUNC) &_RSTr_cpp_rmvnorm, 2},
    {"_RSTr_geteig", (DL_FUNC) &_RSTr_geteig, 1},
    {"_RSTr_m_update_beta", (DL_FUNC) &_RSTr_m_update_beta, 5},
    {"_RSTr_m_update_Z", (DL_FUNC) &_RSTr_m_update_Z, 9},
    {"_RSTr_m_update_G", (DL_FUNC) &_RSTr_m_update_G, 6},
    {"_RSTr_m_update_tau2", (DL_FUNC) &_RSTr_m_update_tau2, 7},
    {"_RSTr_m_update_theta", (DL_FUNC) &_RSTr_m_update_theta, 10},
    {"_RSTr_mst_update_beta", (DL_FUNC) &_RSTr_mst_update_beta, 5},
    {"_RSTr_mst_update_Z", (DL_FUNC) &_RSTr_mst_update_Z, 10},
    {"_RSTr_mst_update_G", (DL_FUNC) &_RSTr_mst_update_G, 7},
    {"_RSTr_mst_update_Ag", (DL_FUNC) &_RSTr_mst_update_Ag, 5},
    {"_RSTr_mst_update_tau2", (DL_FUNC) &_RSTr_mst_update_tau2, 7},
    {"_RSTr_mst_update_theta", (DL_FUNC) &_RSTr_mst_update_theta, 10},
    {"_RSTr_mst_update_rho", (DL_FUNC) &_RSTr_mst_update_rho, 9},
    {"_RSTr_u_update_Z", (DL_FUNC) &_RSTr_u_update_Z, 9},
    {"_RSTr_u_update_sig2", (DL_FUNC) &_RSTr_u_update_sig2, 13},
    {"_RSTr_u_update_tau2", (DL_FUNC) &_RSTr_u_update_tau2, 12},
    {"_RSTr_u_update_theta", (DL_FUNC) &_RSTr_u_update_theta, 10},
    {"_RSTr_u_update_beta", (DL_FUNC) &_RSTr_u_update_beta, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_RSTr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
