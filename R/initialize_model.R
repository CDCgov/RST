#' Initialize CAR model
#' @param name Name of model and corresponding folder
#' @param dir Directory where model will live
#' @param data Dataset including mortality (Y) and population (n) information
#' @param adjacency Dataset including adjacency information
#' @param inits Optional list of initial conditions for each parameter
#' @param priors Optional list of priors for updates
#' @param method Run model with either Binomial data or Poisson data.
#' @param rho_up Controls whether rho update is performed
#' @param impute_lb If counts between lb and ub are suppressed for privacy reasons, impute_lb is lower bound
#' @param impute_ub If counts between lb and ub are suppressed for privacy reasons, impute_ub is upper bound
#' @param seed Set of random seeds to use for data replication
#' @param .ignore_checks If set to TRUE, ignores data checks. Only use if you are certain that your input data is
#' correct and you are encountering bugs during setup.
#'
#' @export
initialize_model = function(
  name,
  dir = getwd(),
  data,
  adjacency,
  inits = NULL,
  priors = NULL,
  method = c("binom", "pois"),
  rho_up = FALSE,
  impute_lb = 1,
  impute_ub = 9,
  seed = 1234,
  .ignore_checks = FALSE
) {
  method = match.arg(method)
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

  num_region = dim(Y)[1]
  num_group = dim(Y)[2]
  num_time = dim(Y)[3]
  set.seed(seed)
  params = list(
    seed = .Random.seed,
    batch = 1,
    total = 0,
    method = method,
    rho_up = rho_up,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    dimnames = dimnames(n)
  )
  spatial_data = get_spatial_data(adjacency, num_region, .ignore_checks)
  priors = get_priors(priors, num_region, num_group, num_time, .ignore_checks)
  inits = get_inits(inits, data, spatial_data$island_id, method, .ignore_checks)
  saveRDS(data, file = paste0(dir, name, "/data.Rds"))
  saveRDS(params, file = paste0(dir, name, "/params.Rds"))
  saveRDS(spatial_data, file = paste0(dir, name, "/spatial_data.Rds"))
  saveRDS(priors, file = paste0(dir, name, "/priors.Rds"))
  saveRDS(inits, file = paste0(dir, name, "/inits.Rds"))
  cat("Model ready!\n")
}
