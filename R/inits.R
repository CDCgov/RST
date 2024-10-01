#' Get initial values
#'
#' @noRd
get_inits = function(inits, data, island_id, method, .ignore_checks) {
  Y = data$Y
  n = data$n
  num_region = dim(Y)[1]
  num_group = dim(Y)[2]
  num_time = dim(Y)[3]
  num_island = length(unique(island_id))
  # Prepare initial values
  initmiss = NULL
  if (is.null(inits$theta)) {
    theta = array(0, dim = c(num_region, num_group, num_time))
    for (time in 1:num_time) {
      for (grp in 1:num_group) {
        if (method == "pois") {
          theta[, grp, time] = log(Y[, grp, time] / n[, grp, time])
          theta[, grp, time] = ifelse(
            !is.finite(theta[, grp, time]),
            log(sum(Y[, grp, time]) / sum(n[, grp, time])),
            theta[, grp, time]
          )
          theta[, grp, time] = ifelse(
            !is.finite(theta[, grp, time]),
            log(sum(Y[, grp, ]) / sum(n[, grp, ])),
            theta[, grp, time]
          )
          theta[, grp, time] = ifelse(
            !is.finite(theta[, grp, time]),
            log(sum(Y, na.rm = TRUE) / sum(n)),
            theta[, grp, time]
          )
        }
        if (method == "binom") {
          theta[, grp, time] = logit(Y[, grp, time] / n[, grp, time])
          theta[, grp, time] = ifelse(
            !is.finite(theta[, grp, time]),
            logit(sum(Y[, grp, time]) / sum(n[, grp, time])),
            theta[, grp, time]
          )
          theta[, grp, time] = ifelse(
            !is.finite(theta[, grp, time]),
            logit(sum(Y[, grp, ]) / sum(n[, grp, ])),
            theta[, grp, time]
          )
          theta[, grp, time] = ifelse(
            !is.finite(theta[, grp, time]),
            logit(sum(Y, na.rm = TRUE) / sum(n)),
            theta[, grp, time]
          )
        }
      }
    }
    inits$theta = theta
    initmiss = c(initmiss, "theta")
  }
  # beta
  if (is.null(inits$beta)) {
    beta = apply(Y, 2:3, sum, na.rm = TRUE) / apply(n, 2:3, sum)
    if (method == "pois") {
      beta = array(log(beta), dim = c(num_island, num_group, num_time))
      beta[!is.finite(beta)] = log(sum(Y, na.rm = TRUE) / sum(n))
    }
    if (method == "binom") {
      beta = array(logit(beta), dim = c(num_island, num_group, num_time))
      beta[!is.finite(beta)] = logit(sum(Y, na.rm = TRUE) / sum(n))
    }
    inits$beta = beta
    initmiss = c(initmiss, "beta")
  }
  # Z
  if (is.null(inits$Z)) {
    inits$Z = inits$theta - inits$beta[island_id + 1, , , drop = FALSE]
    initmiss = c(initmiss, "Z")
  }
  # G
  if (is.null(inits$G)) {
    inits$G = array(diag(num_group) / 7, dim = c(num_group, num_group, num_time))
    initmiss = c(initmiss, "G")
  }
  # rho
  if (is.null(inits$rho)) {
    inits$rho = rep(0.95, num_group)
    initmiss = c(initmiss, "rho")
  }
  # tau2
  if (is.null(inits$tau2)) {
    inits$tau2 = rep(1 / 100, num_group)
    initmiss = c(initmiss, "tau2")
  }
  # Ag
  if (is.null(inits$Ag)) {
    inits$Ag = diag(1 / 7, num_group)
    initmiss = c(initmiss, "Ag")
  }
  if (!.ignore_checks) {
    check_inits(inits, num_region, num_group, num_time, num_island)
  }
  if (!is.null(initmiss)) {
    cat("The following objects were created using defaults in 'inits':", paste(initmiss, collapse = " "), "\n")
  }
  inits
}

#' Check initial values
#'
#' @noRd
#'
check_inits = function(inits, num_region, num_group, num_time, num_island) {
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
