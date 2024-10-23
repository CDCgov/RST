#' Check initial values MSTCAR
#'
#' @noRd
#'
check_inits_u = function(inits, num_region, num_island) {
  cat("Checking inits...\n")
  theta = inits$theta
  beta  = inits$beta
  sig2  = inits$sig2
  tau2  = inits$tau2
  Z     = inits$Z
  chk   = c("theta", "beta", "tau2", "sig2", "Z")
  miss  = sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'inits'
  chk_elem = which(!(names(inits) %in% chk))
  if (length(chk_elem)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'inits':", paste(names(inits)[chk_elem], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'inits':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout = NULL
  errct  = 0
  # theta
  # length doesn't match num_region
  if (length(theta) != num_region) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta is not a length num_region vector. Ensure length(theta) == length(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(theta))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta contains infinite values. Ensure all(is.finite(theta)) or use default value")
    errout = c(errout, errtxt)
  }

  # beta
  # dimensions don't match num_time num_island num_group
  if (length(beta) != num_island) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta is not a length num_island vector. Ensure length(beta) == num_island or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errct  = errct + 1
    errtxt = paste(errct, ": beta contains infinite values. Ensure all(is.finite(beta)) or use default value")
    errout = c(errout, errtxt)
  }
  # sig2
  # is non-positive or infinite
  if ((sig2 <= 0) | !is.finite(sig2)) {
    errct  = errct + 1
    errtxt = paste(errct, ": sig2 is a non-positive or infinite value. Ensure sig2 > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # tau2
  # is non-positive or infinite
  if ((tau2 <= 0) | !is.finite(tau2)) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau2 is a non-positive or infinite value. Ensure tau2 > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # Z
  # dimensions don't match num_region num_group num_time
  if (length(Z) != num_region) {
    errct  = errct + 1
    errtxt = paste(errct, ": Z is not a length num_region vector. Ensure length(Z) == length(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Z contains infinite values. Ensure all(is.finite(Z)) or use default value")
    errout = c(errout, errtxt)
  }
  if (errct) {
    stop(paste(errct, "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
  }
}

#' Check initial values MSTCAR
#'
#' @noRd
#'
check_inits_m = function(inits, num_region, num_group, num_island) {
  cat("Checking inits...\n")
  theta = inits$theta
  beta  = inits$beta
  G     = inits$G
  tau2  = inits$tau2
  Z     = inits$Z
  chk   = c("theta", "beta", "tau2", "G", "Z")
  miss  = sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'inits'
  chk_elem = which(!(names(inits) %in% chk))
  if (length(chk_elem)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'inits':", paste(names(inits)[chk_elem], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'inits':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout = NULL
  errct  = 0
  # theta
  # dimensions don't match num_region num_group
  if (!all(dim(theta) == c(num_region, num_group))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta is not an num_region x num_group array. Ensure dim(theta) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(theta))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta contains infinite values. Ensure all(is.finite(theta)) or use default value")
    errout = c(errout, errtxt)
  }

  # beta
  # dimensions don't match num_time num_island num_group
  if (!all(dim(beta) == c(num_island, num_group))) {
    errct  = errct + 1
    errtxt = paste(errct, ": beta is not an num_island x num_group array. Ensure dim(beta) == num_island x num_group or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errct  = errct + 1
    errtxt = paste(errct, ": beta contains infinite values. Ensure all(is.finite(beta)) or use default value")
    errout = c(errout, errtxt)
  }
  # G
  sig2 = diag(G)
  gcor = G[lower.tri(G)]
  # G not symmetric
  if (!isSymmetric(G)) {
    errct  = errct + 1
    errtxt = paste(errct, ": G is not symmetric. Ensure G is symmetric or use default value")
    errout = c(errout, errtxt)
  }

  # sig2
  # is non-positive or infinite
  if (any((sig2 <= 0) | !is.finite(sig2))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Diagonals of G contain non-positive values. Ensure all diag(G) > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # gcor
  # values are infinite
  if (any(!is.finite(gcor))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Off-diagonals of G contain infinite values. Ensure all(is.finite(G)) or use default value")
    errout = c(errout, errtxt)
  }

  # tau2
  # length not num_group
  if (length(tau2) != num_group) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau2 is not length num_group. Ensure length(tau2) == num_group or use default value")
    errout = c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((tau2 <= 0) | !is.finite(tau2))) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau2 contains non-positive values. Ensure all(tau2 > 0) and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # Z
  # dimensions don't match num_region num_group
  if (!all(dim(Z) == c(num_region, num_group))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Z is not an num_region x num_group array. Ensure dim(Z) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Z contains infinite values. Ensure all(is.finite(Z)) or use default value")
    errout = c(errout, errtxt)
  }

  if (errct) {
    stop(paste(errct, "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
  }
}


#' Check initial values MSTCAR
#'
#' @noRd
#'
check_inits_mst = function(inits, num_region, num_group, num_time, num_island) {
  cat("Checking inits...\n")
  theta = inits$theta
  beta  = inits$beta
  G     = inits$G
  tau2  = inits$tau2
  Z     = inits$Z
  rho   = inits$rho
  Ag    = inits$Ag
  chk   = c("theta", "beta", "tau2", "G", "Ag", "Z", "rho")
  miss  = sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'inits'
  chk_elem = which(!(names(inits) %in% chk))
  if (length(chk_elem)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'inits':", paste(names(inits)[chk_elem], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'inits':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout = NULL
  errct  = 0
  # theta
  # dimensions don't match num_region num_group num_time
  if (!all(dim(theta) == c(num_region, num_group, num_time))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta is not an num_region x num_group x num_time array. Ensure dim(theta) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(theta))) {
    errct  = errct + 1
    errtxt = paste(errct, ": theta contains infinite values. Ensure all(is.finite(theta)) or use default value")
    errout = c(errout, errtxt)
  }

  # beta
  # dimensions don't match num_time num_island num_group
  if (!all(dim(beta) == c(num_island, num_group, num_time))) {
    errct  = errct + 1
    errtxt = paste(errct, ": beta is not an num_island x num_group x num_time array. Ensure dim(beta) == num_island x num_group x num_time or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errct  = errct + 1
    errtxt = paste(errct, ": beta contains infinite values. Ensure all(is.finite(beta)) or use default value")
    errout = c(errout, errtxt)
  }
  # G
  sig2 = apply(G, 3, diag)
  gcor = apply(G, 3, \(G) G[lower.tri(G)])
  # G not symmetric
  #if (any(!sapply(1:num_group, \(j) isSymmetric(G[, , j])))) {
  #  errct  = errct + 1
  #  errtxt = paste(errct, ": Some slices of G are not symmetric. Ensure all slices of G are symmetric or use default value")
  #  errout = c(errout, errtxt)
  #}
  # sig2
  # is non-positive or infinite
  if (any((sig2 <= 0) | !is.finite(sig2))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Diagonals of G contain non-positive values. Ensure all diag(G) > 0 and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # gcor
  # values are infinite
  if (any(!is.finite(gcor))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Off-diagonals of G contain infinite values. Ensure all(is.finite(G)) or use default value")
    errout = c(errout, errtxt)
  }

  # tau2
  # length not num_group
  if (length(tau2) != num_group) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau2 is not length num_group. Ensure length(tau2) == num_group or use default value")
    errout = c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((tau2 <= 0) | !is.finite(tau2))) {
    errct  = errct + 1
    errtxt = paste(errct, ": tau2 contains non-positive values. Ensure all(tau2 > 0) and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # rho
  # length not num_group
  if (length(rho) != num_group) {
    errct  = errct + 1
    errtxt = paste(errct, ": rho is not length num_group. Ensure length(rho) == num_group or use default value")
    errout = c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((rho <= 0) | !is.finite(rho))) {
    errct  = errct + 1
    errtxt = paste(errct, ": rho contains non-positive values. Ensure all(rho > 0) and not infinite or use default value")
    errout = c(errout, errtxt)
  }

  # Z
  # dimensions don't match num_region num_group num_time
  if (!all(dim(Z) == c(num_region, num_group, num_time))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Z is not an num_region x num_group x num_time array. Ensure dim(Z) == dim(Y) or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Z contains infinite values. Ensure all(is.finite(Z)) or use default value")
    errout = c(errout, errtxt)
  }

  # Ag
  # dimensions don't match num_group num_group
  if (!all(dim(Ag) == c(num_group, num_group))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag is not an num_group x num_group matrix. Ensure dim(Ag) == num_group x num_group or use default value")
    errout = c(errout, errtxt)
  }
  # matrix is not symmetric
  if (!isSymmetric(Ag)) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag is not symmetric. Ensure Ag is symmetric or use default value")
    errout = c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Ag))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Ag contains infinite values. Ensure Ag is finite or use default value")
    errout = c(errout, errtxt)
  }
  # diagonals are not positive
  if (any(diag(Ag) <= 0)) {
    errct  = errct + 1
    errtxt = paste(errct, ": diag(Ag) contains non-positive values. Ensure diag(Ag) is positive or use default value")
    errout = c(errout, errtxt)
  }

  if (errct) {
    stop(paste(errct, "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
  }
}
