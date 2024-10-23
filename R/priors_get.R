#' Get priors UCAR
#'
#' @noRd
get_priors_u = function(priors, num_region, .ignore_checks) {
  # Prepare Hyperparameters
  primiss = NULL
  if (is.null(priors$theta_sd)) {
    priors$theta_sd = rep(0.025, num_region)
    primiss = c(primiss, "theta_sd")
  }
  priors$t_accept = rep(0.4, num_region)
  if (is.null(priors$tau_a)) {
    priors$tau_a = 0.001
    primiss = c(primiss, "tau_a")
  }
  if (is.null(priors$tau_b)) {
    priors$tau_b = 0.001
    primiss = c(primiss, "tau_b")
  }
  if (is.null(priors$sig_a)) {
    priors$sig_a = 0.001
    primiss = c(primiss, "sig_a")
  }
  if (is.null(priors$sig_b)) {
    priors$sig_b = 0.001
    primiss = c(primiss, "sig_b")
  }
  if (!.ignore_checks) {
    check_priors_u(priors, num_region)
  }
  if (!is.null(primiss)) {
    cat("The following objects were created using defaults in 'priors':", paste(primiss, collapse = " "), "\n")
  }
  priors
}

#' Get priors MCAR
#'
#' @noRd
get_priors_m = function(priors, num_region, num_group, .ignore_checks) {
  # Prepare Hyperparameters
  primiss = NULL
  if (is.null(priors$theta_sd)) {
    priors$theta_sd = array(0.025, dim = c(num_region, num_group))
    primiss = c(primiss, "theta_sd")
  }
  priors$t_accept = array(0.4, dim = c(num_region, num_group))
  if (is.null(priors$tau_a)) {
    priors$tau_a = 0.001
    primiss = c(primiss, "tau_a")
  }
  if (is.null(priors$tau_b)) {
    priors$tau_b = 0.001
    primiss = c(primiss, "tau_b")
  }
  if (is.null(priors$G_scale)) {
    priors$G_scale = diag(1 / 7, num_group)
    primiss = c(primiss, "G_scale")
  }
  if (is.null(priors$G_df)) {
    priors$G_df = num_group + 2
    primiss = c(primiss, "G_df")
  }
  if (!.ignore_checks) {
    check_priors_m(priors, num_region, num_group)
  }
  if (!is.null(primiss)) {
    cat("The following objects were created using defaults in 'priors':", paste(primiss, collapse = " "), "\n")
  }
  priors
}

#' Get priors USTCAR
#'
#' @noRd
get_priors_ust = function(priors, num_region, num_time, .ignore_checks) {
  primiss = NULL
  if (is.null(priors$theta_sd)) {
    priors$theta_sd = array(0.025, dim = c(num_region, num_time))
    primiss = c(primiss, "theta_sd")
  }
  priors$t_accept = array(0.4, dim = c(num_region, num_time))
  #if (!.ignore_checks) {
  #  check_priors_ust(priors, num_region, num_time)
  #}
  if (!is.null(primiss)) {
    cat("The following objects were created using defaults in 'priors':", paste(primiss, collapse = " "), "\n")
  }
  priors
}

#' Get priors MSTCAR
#'
#' @noRd
get_priors_mst = function(priors, num_region, num_group, num_time, .ignore_checks) {
  # Prepare Hyperparameters
  primiss = NULL
  if (is.null(priors$theta_sd)) {
    priors$theta_sd = array(0.025, dim = c(num_region, num_group, num_time))
    primiss = c(primiss, "theta_sd")
  }
  priors$t_accept = array(0.4, dim = c(num_region, num_group, num_time))
  if (is.null(priors$tau_a)) {
    priors$tau_a = 0.001
    primiss = c(primiss, "tau_a")
  }
  if (is.null(priors$tau_b)) {
    priors$tau_b = 0.001
    primiss = c(primiss, "tau_b")
  }
  if (is.null(priors$Ag_scale)) {
    priors$Ag_scale = diag(1 / 7, num_group)
    primiss = c(primiss, "Ag_scale")
  }
  if (is.null(priors$Ag_df)) {
    priors$Ag_df = num_group + 2
    primiss = c(primiss, "Ag_df")
  }
  if (is.null(priors$G_df)) {
    priors$G_df = num_group + 2
    primiss = c(primiss, "G_df")
  }
  if (is.null(priors$rho_a)) {
    priors$rho_a = 95
    primiss = c(primiss, "rho_a")
  }
  if (is.null(priors$rho_b)) {
    priors$rho_b = 5
    primiss = c(primiss, "rho_b")
  }
  if (is.null(priors$rho_sd)) {
    priors$rho_sd = rep(0.05, num_group)
    primiss = c(primiss, "rho_sd")
  }
  priors$r_accept = rep(0.4, num_group)
  if (!.ignore_checks) {
    check_priors_mst(priors, num_region, num_group, num_time)
  }
  if (!is.null(primiss)) {
    cat("The following objects were created using defaults in 'priors':", paste(primiss, collapse = " "), "\n")
  }
  priors
}
