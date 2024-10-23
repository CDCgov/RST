#' Check priors MSTCAR
#'
#' @noRd
check_priors_u = function(priors, num_region) {
  cat("Checking priors...\n")
  tau_a    = priors$tau_a
  tau_b    = priors$tau_b
  sig_a    = priors$sig_a
  sig_b    = priors$sig_b
  theta_sd = priors$theta_sd
  chk  = c("tau_a", "tau_b", "sig_a", "sig_b", "theta_sd")
  miss = sapply(1:length(chk), \(x) !any(names(priors) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'priors': ", paste(chk[miss], collapse = ", "))
  }
  chk = c(chk, "t_accept")
  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'priors'
  chk_elem = which(!(names(priors) %in% chk))
  if (length(chk_elem)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'priors':", paste(names(priors)[chk_elem], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'priors':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout = NULL
  errct  = 0

  # tau_a
  # is non-positive or infinite
  if ((tau_a <= 0) | !is.finite(tau_a)) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau_a is not positive. Ensure tau_a > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # tau_b
  # is non-positive or infinite
  if ((tau_b <= 0) | !is.finite(tau_b)) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau_b is not positive. Ensure tau_b > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # sig_a
  # is non-positive or infinite
  if ((sig_a <= 0) | !is.finite(sig_a)) {
    errct  = errct + 1
    errtxt = paste(errct, ": sig_a is not positive. Ensure sig_a > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # sig_b
  # is non-positive or infinite
  if ((sig_b <= 0) | !is.finite(sig_b)) {
    errct  = errct + 1
    errtxt = paste(errct, ": sig_b is not positive. Ensure sig_b > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # theta_sd
  # dim not num_time num_region num_group
  if (length(theta_sd) != num_region) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta_sd has different length than data. Ensure length(theta_sd) == length(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((theta_sd <= 0) | !is.finite(theta_sd))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta_sd contains non-positive values. Ensure all(theta_sd > 0) and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  if (errct) {
    stop(paste(errct, "error(s) found in list 'priors':\n", paste(errout, collapse = "\n ")))
  }
}


#' Check priors MCAR
#'
#' @noRd
check_priors_m = function(priors, num_region, num_group, num_time) {
  cat("Checking priors...\n")
  G_scale  = priors$G_scale
  G_df     = priors$G_df
  tau_a    = priors$tau_a
  tau_b    = priors$tau_b
  theta_sd = priors$theta_sd
  chk  = c("G_scale", "G_df", "tau_a", "tau_b", "theta_sd")
  miss = sapply(1:length(chk), \(x) !any(names(priors) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'priors': ", paste(chk[miss], collapse = ", "))
  }
  chk = c(chk, "t_accept")
  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'priors'
  chk_elem = which(!(names(priors) %in% chk))
  if (length(chk_elem)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'priors':", paste(names(priors)[chk_elem], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'priors':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout = NULL
  errct  = 0
  # G_scale
  # dimensions don't match num_group num_group
  if (!all(dim(G_scale) == c(num_group, num_group))) {
    errct  = errct + 1
    errtxt = paste(errct, ": G_scale is not an num_group x num_group matrix. Ensure dim(G_scale) == num_group x num_group or use default value")
    errout = c(errout, errtxt)
  }
  # matrix is not symmetric
  if (!isSymmetric(G_scale)) {
    errct  = errct + 1
    errtxt = paste(errct, ": G_scale is not symmetric. Ensure G_scale is symmetric or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(G_scale))) {
    errct  = errct + 1
    errtxt = paste(errct, ": G_scale contains infinite values. Ensure G_scale is finite or use default value")
    errout = c(errout, errtxt)
  }
  # diagonals are not positive
  if (any(diag(G_scale) <= 0)) {
    errct  = errct + 1
    errtxt = paste(errct, ": diag(G_scale) contains non-positive values. Ensure diag(G_scale) is positive or use default value")
    errout = c(errout, errtxt)
  }

  # G_df
  # is not a whole number
  if (floor(G_df) != G_df) {
    errct  = errct + 1
    errtxt = paste(errct, ": G_df is not a whole number. Ensure G_df is whole number or use default value")
    errout = c(errout, errtxt)
  }
  # is less than df
  if (G_df <= num_group - 1) {
    errct  = errct + 1
    errtxt = paste(errct, ": G_df too small. Ensure G_df > num_group - 1 or use default value")
    errout = c(errout, errtxt)
  }

  # tau_a
  # is non-positive or infinite
  if ((tau_a <= 0) | !is.finite(tau_a)) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau_a is not positive. Ensure tau_a > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # tau_b
  # is non-positive or infinite
  if ((tau_b <= 0) | !is.finite(tau_b)) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau_b is not positive. Ensure tau_b > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # theta_sd
  # dim not num_region num_group
  if (!all(dim(theta_sd) == c(num_region, num_group))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta_sd has different dimensions than data. Ensure dim(theta_sd) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((theta_sd <= 0) | !is.finite(theta_sd))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta_sd contains non-positive values. Ensure all(theta_sd > 0) and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  if (errct) {
    stop(paste(errct, "error(s) found in list 'priors':\n", paste(errout, collapse = "\n ")))
  }
}



#' Check priors MSTCAR
#'
#' @noRd
check_priors_mst = function(priors, num_region, num_group, num_time) {
  cat("Checking priors...\n")
  Ag_scale = priors$Ag_scale
  Ag_df    = priors$Ag_df
  G_df     = priors$G_df
  tau_a    = priors$tau_a
  tau_b    = priors$tau_b
  rho_a    = priors$rho_a
  rho_b    = priors$rho_b
  theta_sd = priors$theta_sd
  rho_sd   = priors$rho_sd
  chk  = c("Ag_scale", "G_df", "Ag_df", "tau_a", "tau_b", "rho_a", "rho_b", "theta_sd", "rho_sd")
  miss = sapply(1:length(chk), \(x) !any(names(priors) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'priors': ", paste(chk[miss], collapse = ", "))
  }
  chk = c(chk, "t_accept", "r_accept")
  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'priors'
  chk_elem = which(!(names(priors) %in% chk))
  if (length(chk_elem)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'priors':", paste(names(priors)[chk_elem], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'priors':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout = NULL
  errct  = 0
  # Ag_scale
  # dimensions don't match num_group num_group
  if (!all(dim(Ag_scale) == c(num_group, num_group))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag_scale is not an num_group x num_group matrix. Ensure dim(Ag_scale) == num_group x num_group or use default value")
    errout = c(errout, errtxt)
  }
  # matrix is not symmetric
  if (!isSymmetric(Ag_scale)) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag_scale is not symmetric. Ensure Ag_scale is symmetric or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Ag_scale))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag_scale contains infinite values. Ensure Ag_scale is finite or use default value")
    errout = c(errout, errtxt)
  }
  # diagonals are not positive
  if (any(diag(Ag_scale) <= 0)) {
    errct  = errct + 1
    errtxt = paste(errct, ": diag(Ag_scale) contains non-positive values. Ensure diag(Ag_scale) is positive or use default value")
    errout = c(errout, errtxt)
  }

  # G_df
  # is not a whole number
  if (floor(G_df) != G_df) {
    errct  = errct + 1
    errtxt = paste(errct, ": G_df is not a whole number. Ensure G_df is whole number or use default value")
    errout = c(errout, errtxt)
  }
  # is less than df
  if (G_df <= num_group - 1) {
    errct  = errct + 1
    errtxt = paste(errct, ": G_df too small. Ensure G_df > num_group - 1 or use default value")
    errout = c(errout, errtxt)
  }

  # Ag_df
  # is not a whole number
  if (floor(Ag_df) != Ag_df) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag_df is not a whole number. Ensure Ag_df is whole number or use default value")
    errout = c(errout, errtxt)
  }
  # is less than df
  if (Ag_df <= num_group - 1) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag_df too small. Ensure Ag_df > num_group - 1 or use default value")
    errout = c(errout, errtxt)
  }

  # tau_a
  # is non-positive or infinite
  if ((tau_a <= 0) | !is.finite(tau_a)) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau_a is not positive. Ensure tau_a > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # tau_b
  # is non-positive or infinite
  if ((tau_b <= 0) | !is.finite(tau_b)) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau_b is not positive. Ensure tau_b > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # rho_a
  # is non-positive or infinite
  if ((rho_a <= 0) | !is.finite(rho_a)) {
    errct  = errct + 1
    errtxt = paste(errct, ": rho_a is not positive. Ensure rho_a > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # rho_b
  # is non-positive or infinite
  if ((rho_b <= 0) | !is.finite(rho_b)) {
    errct  = errct + 1
    errtxt = paste(errct, ": rho_b is not positive. Ensure rho_b > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # theta_sd
  # dim not num_time num_region num_group
  if (!all(dim(theta_sd) == c(num_region, num_group, num_time))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta_sd has different dimensions than data. Ensure dim(theta_sd) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((theta_sd <= 0) | !is.finite(theta_sd))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta_sd contains non-positive values. Ensure all(theta_sd > 0) and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # rho_sd
  # length not num_group
  if (length(rho_sd) != num_group) {
    errct  = errct + 1
    errtxt = paste(errct, ": rho_sd is not length num_group. Ensure length(rho_sd) == num_group or use default value")
    errout = c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((rho_sd <= 0) | !is.finite(rho_sd))) {
    errct  = errct + 1
    errtxt = paste(errct, ": rho_sd contains non-positive values. Ensure all(rho_sd > 0) and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  if (errct) {
    stop(paste(errct, "error(s) found in list 'priors':\n", paste(errout, collapse = "\n ")))
  }
}
