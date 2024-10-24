#' Initialize CAR model
#' @param name Name of model and corresponding folder
#' @param dir Directory where model will live
#' @param data Dataset including mortality (Y) and population (n) information
#' @param adjacency Dataset including adjacency information
#' @param inits Optional list of initial conditions for each parameter
#' @param priors Optional list of priors for updates
#' @param model Run model as an MSTCAR/UCAR/USTCAR/MCAR model
#' @param method Run model with either Binomial data or Poisson data
#' @param m0 Baseline neighbor count by region
#' @param A Describes intensity of smoothing between regions
#' @param rho_up Controls whether rho update is performed for USTCAR or MSTCAR models
#' @param impute_lb If counts between lb and ub are suppressed for privacy reasons, impute_lb is lower bound
#' @param impute_ub If counts between lb and ub are suppressed for privacy reasons, impute_ub is upper bound
#' @param seed Set of random seeds to use for data replication
#' @param .ignore_checks If set to TRUE, ignores data checks. Only use if you are certain that your input data is
#' correct and you are encountering bugs during setup
#'
#' @export
initialize_model = function(
  name,
  dir = getwd(),
  data,
  adjacency,
  inits = NULL,
  priors = NULL,
  model = c("mstcar", "ucar", "ustcar", "mcar"),
  method = c("binom", "pois"),
  m0 = 3,
  A = NULL,
  rho_up = FALSE,
  impute_lb = 1,
  impute_ub = 9,
  seed = 1234,
  .ignore_checks = FALSE
) {
  method = match.arg(method)
  model = match.arg(model)
  if (model == "ucar") {
    initialize_model_u(name, dir, data, adjacency, inits, priors, method, m0, A, impute_lb, impute_ub, seed, .ignore_checks)
  }
  if (model == "mcar") {
    initialize_model_m(name, dir, data, adjacency, inits, priors, method, m0, A, impute_lb, impute_ub, seed, .ignore_checks)
  }
  #if (model == "ustcar") {
  #  initialize_model_ust(name, dir, data, adjacency, inits, priors, method, m0, A, rho_up, impute_lb, impute_ub, seed, .ignore_checks)
  #}
  if (model == "mstcar") {
    initialize_model_mst(name, dir, data, adjacency, inits, priors, method, m0, A, rho_up, impute_lb, impute_ub, seed, .ignore_checks)
  }
}

#' Initialize UCAR model
#'
#' @noRd
initialize_model_u = function(name, dir, data, adjacency, inits, priors, method, m0, A, impute_lb, impute_ub, seed, .ignore_checks) {
  Y = data$Y
  n = data$n
  if (!.ignore_checks) {
    check_data(data)
  }
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  if (!dir.exists(paste0(dir, name))) {
    dir.create(paste0(dir, name))
  }
  pars = c("theta", "beta", "Z", "sig2", "tau2")
  for (par in pars) {
    if (!dir.exists(paste0(dir, name, "/", par))) {
      dir.create(paste0(dir, name, "/", par))
    }
  }
  if (is.null(A)) {
    A = 1
  }

  num_region = length(Y)
  set.seed(seed)
  params = list(
    seed = .Random.seed,
    batch = 1,
    total = 0,
    model = "ucar",
    method = method,
    m0 = m0,
    A = A,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    dimnames = names(n)
  )
  spatial_data = get_spatial_data(adjacency, num_region, .ignore_checks)
  priors = get_priors_u(priors, num_region, .ignore_checks)
  inits = get_inits_u(inits, data, spatial_data$island_id, method, .ignore_checks)
  saveRDS(data, file = paste0(dir, name, "/data.Rds"))
  saveRDS(params, file = paste0(dir, name, "/params.Rds"))
  saveRDS(spatial_data, file = paste0(dir, name, "/spatial_data.Rds"))
  saveRDS(priors, file = paste0(dir, name, "/priors.Rds"))
  saveRDS(inits, file = paste0(dir, name, "/inits.Rds"))
  cat("Model ready!\n")
}

#' Initialize MCAR model
#'
#' @noRd
initialize_model_m = function(name, dir, data, adjacency, inits, priors, method, m0, A, impute_lb, impute_ub, seed, .ignore_checks) {
  Y = data$Y
  n = data$n
  if (!.ignore_checks) {
    check_data(data)
  }
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  if (!dir.exists(paste0(dir, name))) {
    dir.create(paste0(dir, name))
  }
  pars = c("theta", "beta", "Z", "G", "tau2")
  for (par in pars) {
    if (!dir.exists(paste0(dir, name, "/", par))) {
      dir.create(paste0(dir, name, "/", par))
    }
  }

  num_region = dim(Y)[1]
  num_group = dim(Y)[2]
  if (is.null(A)) {
    A = num_group * colSums(Y) / sum(Y)
  }
  set.seed(seed)
  params = list(
    seed = .Random.seed,
    batch = 1,
    total = 0,
    model = "mcar",
    method = method,
    m0 = m0,
    A = A,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    dimnames = dimnames(n)
  )
  spatial_data = get_spatial_data(adjacency, num_region, .ignore_checks)
  priors = get_priors_m(priors, num_region, num_group, .ignore_checks)
  inits = get_inits_m(inits, data, spatial_data$island_id, method, .ignore_checks)
  saveRDS(data, file = paste0(dir, name, "/data.Rds"))
  saveRDS(params, file = paste0(dir, name, "/params.Rds"))
  saveRDS(spatial_data, file = paste0(dir, name, "/spatial_data.Rds"))
  saveRDS(priors, file = paste0(dir, name, "/priors.Rds"))
  saveRDS(inits, file = paste0(dir, name, "/inits.Rds"))
  cat("Model ready!\n")
}

#' Initialize USTCAR model
#'
#' @noRd
initialize_model_ust = function(name, dir, data, adjacency, inits, priors, method, m0, A, rho_up, impute_lb, impute_ub, seed, .ignore_checks) {
  Y = data$Y
  n = data$n
  if (!.ignore_checks) {
    check_data(data)
  }
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  if (!dir.exists(paste0(dir, name))) {
    dir.create(paste0(dir, name))
  }
  pars = c("theta", "beta", "Z", "sig2", "tau2")
  if (rho_up) {
    pars = c(pars, "rho")
  }
  for (par in pars) {
    if (!dir.exists(paste0(dir, name, "/", par))) {
      dir.create(paste0(dir, name, "/", par))
    }
  }
  # Graph of total cases
  plot(dimnames(Y)[[2]], apply(Y, 2, sum, na.rm = TRUE))
  plot(dimnames(Y)[[2]], apply(n, 2, sum))


  num_region = dim(Y)[1]
  num_time = dim(Y)[2]
  set.seed(seed)
  params = list(
    seed = .Random.seed,
    batch = 1,
    total = 0,
    model = "ustcar",
    method = method,
    rho_up = rho_up,
    m0 = m0,
    A = A,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    dimnames = dimnames(n)
  )
  spatial_data = get_spatial_data(adjacency, num_region, .ignore_checks)
  priors = get_priors_ust(priors, num_region, num_time, .ignore_checks)
  inits = get_inits_ust(inits, data, spatial_data$island_id, method, .ignore_checks)
  saveRDS(data, file = paste0(dir, name, "/data.Rds"))
  saveRDS(params, file = paste0(dir, name, "/params.Rds"))
  saveRDS(spatial_data, file = paste0(dir, name, "/spatial_data.Rds"))
  saveRDS(priors, file = paste0(dir, name, "/priors.Rds"))
  saveRDS(inits, file = paste0(dir, name, "/inits.Rds"))
  cat("Model ready!\n")
}

#' Initialize MSTCAR model
#'
#' @noRd
initialize_model_mst = function(name, dir, data, adjacency, inits, priors, method, m0, A, rho_up, impute_lb, impute_ub, seed, .ignore_checks) {
  Y = data$Y
  n = data$n
  if (!.ignore_checks) {
    check_data(data)
  }
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  if (!dir.exists(paste0(dir, name))) {
    dir.create(paste0(dir, name))
  }
  pars = c("theta", "beta", "Z", "G", "Ag", "tau2")
  if (rho_up) {
    pars = c(pars, "rho")
  }
  for (par in pars) {
    if (!dir.exists(paste0(dir, name, "/", par))) {
      dir.create(paste0(dir, name, "/", par))
    }
  }
  # Graph of total cases
  plot(dimnames(Y)[[3]], apply(Y, 3, sum, na.rm = TRUE))
  plot(dimnames(Y)[[3]], apply(n, 3, sum))

  if (is.null(A)) {
    A = NULL
  }
  num_region = dim(Y)[1]
  num_group = dim(Y)[2]
  num_time = dim(Y)[3]
  set.seed(seed)
  params = list(
    seed = .Random.seed,
    batch = 1,
    total = 0,
    model = "mstcar",
    method = method,
    rho_up = rho_up,
    m0 = m0,
    A = A,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    dimnames = dimnames(n)
  )
  spatial_data = get_spatial_data(adjacency, num_region, .ignore_checks)
  priors = get_priors_mst(priors, num_region, num_group, num_time, .ignore_checks)
  inits = get_inits_mst(inits, data, spatial_data$island_id, method, .ignore_checks)
  saveRDS(data, file = paste0(dir, name, "/data.Rds"))
  saveRDS(params, file = paste0(dir, name, "/params.Rds"))
  saveRDS(spatial_data, file = paste0(dir, name, "/spatial_data.Rds"))
  saveRDS(priors, file = paste0(dir, name, "/priors.Rds"))
  saveRDS(inits, file = paste0(dir, name, "/inits.Rds"))
  cat("Model ready!\n")
}
