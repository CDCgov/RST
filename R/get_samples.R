#' Load MCMC samples
#' @param name  Name of model
#' @param dir   Directory where model lives
#' @param param Which parameter samples to load
#' @param burn  Numer of burn-in samples to discard
#'
#' @export
load_samples = function(
  name,
  dir,
  param = c("theta", "beta", "Z", "G", "Ag", "tau2", "rho"),
  burn = 2000
) {
  param  = match.arg(param)
  mar    = c("theta" = 4, "beta" = 4, "Z" = 4, "G" = 4, "Ag" = 3, "tau2" = 2, "rho" = 2)
  #batch  = which(c(1:20 * 50, 1:20 * 250 + 1000) > burn)
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
      dimnames(output) = list(island = as.character(1:num_island), group = dims[[2]], time = dims[[3]], its = its)
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
