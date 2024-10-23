#' Get initial values UCAR
#'
#' @noRd
get_inits_u = function(inits, data, island_id, method, .ignore_checks) {
  Y = data$Y
  n = data$n
  num_region = length(Y)
  num_island = length(unique(island_id))
  # Prepare initial values
  initmiss = NULL
  if (is.null(inits$theta)) {
    if (method == "pois") {
      theta = log(Y / n)
      theta[!is.finite(theta)] = log(sum(Y, na.rm = TRUE) / sum(n))
    }
    if (method == "binom") {
      theta = logit(Y / n)
      theta[!is.finite(theta)] = logit(sum(Y, na.rm = TRUE) / sum(n))
    }
    inits$theta = theta
    initmiss = c(initmiss, "theta")
  }
  # beta
  if (is.null(inits$beta)) {
    beta = sum(Y, na.rm = TRUE) / sum(n)
    if (method == "pois") {
      beta = rep(log(beta), num_island)
    }
    if (method == "binom") {
      beta = rep(logit(beta), num_island)
    }
    inits$beta = beta
    initmiss = c(initmiss, "beta")
  }
  # Z
  if (is.null(inits$Z)) {
    inits$Z = inits$theta - inits$beta[island_id + 1]
    initmiss = c(initmiss, "Z")
  }
  # tau2
  if (is.null(inits$tau2)) {
    inits$tau2 = 1 / 100
    initmiss = c(initmiss, "tau2")
  }
  # sig2
  if (is.null(inits$sig2)) {
    inits$sig2 = 1 / 100
    initmiss = c(initmiss, "sig2")
  }
  if (!.ignore_checks) {
    check_inits_u(inits, num_region, num_island)
  }
  if (!is.null(initmiss)) {
    cat("The following objects were created using defaults in 'inits':", paste(initmiss, collapse = " "), "\n")
  }
  inits
}

#' Get initial values MCAR
#'
#' @noRd
get_inits_m = function(inits, data, island_id, method, .ignore_checks) {
  Y = data$Y
  n = data$n
  num_region = dim(Y)[1]
  num_group = dim(Y)[2]
  num_island = length(unique(island_id))
  # Prepare initial values
  initmiss = NULL
  # beta
  if (is.null(inits$beta)) {
    beta = apply(Y, 2, sum, na.rm = TRUE) / apply(n, 2, sum)
    if (method == "pois") {
      beta = t(array(log(beta), dim = c(num_group, num_island)))
      beta[!is.finite(beta)] = log(sum(Y, na.rm = TRUE) / sum(n))
    }
    if (method == "binom") {
      beta = t(array(logit(beta), dim = c(num_group, num_island)))
      beta[!is.finite(beta)] = logit(sum(Y, na.rm = TRUE) / sum(n))
    }
    inits$beta = beta
    initmiss = c(initmiss, "beta")
  }
  if (is.null(inits$theta)) {
    if (method == "pois") {
      theta = log(Y / n)
      theta[!is.finite(theta)] = beta[island_id + 1, ][which(!is.finite(theta))]
    }
    if (method == "binom") {
      theta = logit(Y / n)
      theta[!is.finite(theta)] = beta[island_id + 1, ][which(!is.finite(theta))]
    }
    inits$theta = theta
    initmiss = c(initmiss, "theta")
  }
  # Z
  if (is.null(inits$Z)) {
    inits$Z = inits$theta - inits$beta[island_id + 1, , drop = FALSE]
    initmiss = c(initmiss, "Z")
  }
  # G
  if (is.null(inits$G)) {
    inits$G = diag(num_group) / 7
    initmiss = c(initmiss, "G")
  }
  # tau2
  if (is.null(inits$tau2)) {
    inits$tau2 = rep(1 / 100, num_group)
    initmiss = c(initmiss, "tau2")
  }
  if (!.ignore_checks) {
    check_inits_m(inits, num_region, num_group, num_island)
  }
  if (!is.null(initmiss)) {
    cat("The following objects were created using defaults in 'inits':", paste(initmiss, collapse = " "), "\n")
  }
  inits
}

#' Get initial values USTCAR
#'
#' @noRd
get_inits_ust = function(inits, data, island_id, method, .ignore_checks) {
  Y = data$Y
  n = data$n
  num_region = dim(Y)[1]
  num_time = dim(Y)[2]
  num_island = length(unique(island_id))
  # Prepare initial values
  initmiss = NULL
  if (is.null(inits$theta)) {
    if (method == "pois") {
      theta = log(Y / n)
      theta[!is.finite(theta)] = log(sum(Y, na.rm = TRUE) / sum(n))
    }
    if (method == "binom") {
      theta = logit(Y / n)
      theta[!is.finite(theta)] = logit(sum(Y, na.rm = TRUE) / sum(n))
    }
    inits$theta = theta
    initmiss = c(initmiss, "theta")
  }
  # beta
  if (is.null(inits$beta)) {
    beta = apply(Y, 2, sum, na.rm = TRUE) / apply(n, 2, sum)
    if (method == "pois") {
      beta = array(log(beta), dim = c(num_island, num_time))
      beta[!is.finite(beta)] = log(sum(Y, na.rm = TRUE) / sum(n))
    }
    if (method == "binom") {
      beta = array(logit(beta), dim = c(num_island, num_time))
      beta[!is.finite(beta)] = logit(sum(Y, na.rm = TRUE) / sum(n))
    }
    inits$beta = beta
    initmiss = c(initmiss, "beta")
  }
  # Z
  if (is.null(inits$Z)) {
    inits$Z = inits$theta - inits$beta[island_id + 1, , drop = FALSE]
    initmiss = c(initmiss, "Z")
  }
  # rho
  if (is.null(inits$rho)) {
    inits$rho = 0.95
    initmiss = c(initmiss, "rho")
  }
  # tau2
  if (is.null(inits$tau2)) {
    inits$tau2 = 1 / 100
    initmiss = c(initmiss, "tau2")
  }
  # sig2?
  if (is.null(inits$sig2)) {

  }
  #if (!.ignore_checks) {
  #  check_inits_ust(inits, num_region, num_time, num_island)
  #}
  if (!is.null(initmiss)) {
    cat("The following objects were created using defaults in 'inits':", paste(initmiss, collapse = " "), "\n")
  }
  inits
}

#' Get initial values MSTCAR
#'
#' @noRd
get_inits_mst = function(inits, data, island_id, method, .ignore_checks) {
  Y = data$Y
  n = data$n
  num_region = dim(Y)[1]
  num_group = dim(Y)[2]
  num_time = dim(Y)[3]
  num_island = length(unique(island_id))
  # Prepare initial values
  initmiss = NULL
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
  # theta
  if (is.null(inits$theta)) {
    if (method == "pois") {
      theta = log(Y / n)
      theta[!is.finite(theta)] = beta[island_id + 1, , ][which(!is.finite(theta))]
    }
    if (method == "binom") {
      theta = logit(Y / n)
      theta[!is.finite(theta)] = beta[island_id + 1, , ][which(!is.finite(theta))]
    }
    inits$theta = theta
    initmiss = c(initmiss, "theta")
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
    check_inits_mst(inits, num_region, num_group, num_time, num_island)
  }
  if (!is.null(initmiss)) {
    cat("The following objects were created using defaults in 'inits':", paste(initmiss, collapse = " "), "\n")
  }
  inits
}
