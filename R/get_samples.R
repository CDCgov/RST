#' Load MCMC samples
#' @param name  Name of model
#' @param dir   Directory where model lives
#' @param time  TIme period of data to pull from. `last` gathers final time period data and `all` gathers all years of data.
#' @param param Which parameter samples to load
#' @param burn  Numer of burn-in samples to discard
#'
#' @export
load_samples = function(name, dir, param = "theta", burn = 2000) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  model = readRDS(paste0(dir, name, "/params.Rds"))$model
  samples = NULL
  if (model == "ucar") {
    samples = load_samples_u(name, dir, param, burn)
  }
  if (model == "mcar") {
    samples = load_samples_m(name, dir, param, burn)
  }
  #if (model == "ustcar") {
  #  samples = load_samples_ust(name, dir, param, burn)
  #}
  if (model == "mstcar") {
    samples = load_samples_mst(name, dir, param, burn)
  }
  samples
}

#' Load MCMC samples, UCAR
#'
#' @noRd
load_samples_u = function(name, dir, param, burn) {
  mar   = c("theta" = 2, "beta" = 2, "Z" = 2, "sig2" = 1, "tau2" = 1)
  batch = which(1:60 * 100 > burn)
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  params = readRDS(paste0(dir, name, "/params.Rds"))
  files  = paste0(dir, name, "/", param, "/", param, "_out_", batch, ".Rds")
  output = abind::abind(lapply(files, readRDS), along = mar[param])
  if (param %in% c("theta", "beta")) {
    if (params$method == "binom") output = expit(output)
    if (params$method == "pois" ) output = exp(output)
  }
  dims = params$dimnames
  if (!is.null(dims)) {
    its = seq(burn + 10, 6000, by = 10)
    if (param == "beta") {
      num_island = readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
      dimnames(output) = list(island = 1:num_island, its = its)
    }
    if (param %in% c("Z", "theta")) {
      dimnames(output) = list(region = dims, its = its)
    }
  }
  output
}

#' Load MCMC samples, MCAR
#'
#' @noRd
load_samples_m = function(name, dir, param, burn) {
  mar   = c("theta" = 3, "beta" = 3, "Z" = 3, "G" = 3, "tau2" = 2)
  batch = which(1:60 * 100 > burn)
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  params = readRDS(paste0(dir, name, "/params.Rds"))
  files  = paste0(dir, name, "/", param, "/", param, "_out_", batch, ".Rds")
  output = abind::abind(lapply(files, readRDS), along = mar[param])
  if (param %in% c("theta", "beta")) {
    if (params$method == "binom") output = expit(output)
    if (params$method == "pois" ) output = exp(output)
  }
  dims = params$dimnames
  if (!is.null(dims)) {
    its = seq(burn + 10, 6000, by = 10)
    if (param == "beta") {
      num_island = readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
      dimnames(output) = list(island = 1:num_island, group = dims[[2]], time = dims[[3]], its = its)
    }
    if (param %in% c("Z", "theta")) {
      dimnames(output) = c(dims, list(its = its))
    }
    if (param %in% c("tau2")) {
      dimnames(output) = list(group = dims[[2]], its = its)
    }
    if (param == "G") {
      dimnames(output) = list(group1 = dims[[2]], group2 = dims[[2]], its = its)
    }
  }
  output
}

#' Load MCMC samples, MSTCAR
#'
#' @noRd
load_samples_mst = function(name, dir, param, burn) {
  mar   = c("theta" = 4, "beta" = 4, "Z" = 4, "G" = 4, "Ag" = 3, "tau2" = 2, "rho" = 2)
  batch = which(1:60 * 100 > burn)
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  params = readRDS(paste0(dir, name, "/params.Rds"))
  files  = paste0(dir, name, "/", param, "/", param, "_out_", batch, ".Rds")
  output = abind::abind(lapply(files, readRDS), along = mar[param])
  if (param %in% c("theta", "beta")) {
    if (params$method == "binom") output = expit(output)
    if (params$method == "pois" ) output = exp(output)
  }
  dims = params$dimnames
  if (!is.null(dims)) {
    its = seq(burn + 10, 6000, by = 10)
    if (param == "beta") {
      num_island = readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
      dimnames(output) = list(island = 1:num_island, group = dims[[2]], time = dims[[3]], its = its)
    }
    if (param %in% c("Z", "theta")) {
      dimnames(output) = c(dims, list(its = its))
    }
    if (param %in% c("tau2", "rho")) {
      dimnames(output) = list(group = dims[[2]], its = its)
    }
    if (param == "Ag") {
      dimnames(output) = list(group1 = dims[[2]], group2 = dims[[2]], its = its)
    }
    if (param == "G") {
      dimnames(output) = list(group1 = dims[[2]], group2 = dims[[2]], time = dims[[3]], its = its)
    }
  }
  output
}
