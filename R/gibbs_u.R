#' Gibbs sampler
#' @useDynLib RST, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#'
#' @noRd
gibbs_u = function(name, dir, .show_plots) {
  data = readRDS(paste0(dir, name, "/data.Rds"))
  Y    = data$Y
  n    = data$n
  miss = which(!is.finite(Y))

  spatial_data = readRDS(paste0(dir, name, "/spatial_data.Rds"))
  adjacency     = spatial_data$adjacency
  num_adj       = spatial_data$num_adj
  island_region = spatial_data$island_region
  island_id     = spatial_data$island_id
  num_island    = spatial_data$num_island
  num_island_region = sapply(island_region, length)

  priors = readRDS(paste0(dir, name, "/priors.Rds"))
  tau_a    = priors$tau_a
  tau_b    = priors$tau_b
  sig_a    = priors$sig_a
  sig_b    = priors$sig_b
  theta_sd = priors$theta_sd
  t_accept = priors$t_accept

  inits = readRDS(paste0(dir, name, "/inits.Rds"))
  theta = inits$theta
  beta  = inits$beta
  Z     = inits$Z
  sig2  = inits$sig2
  tau2  = inits$tau2

  params = readRDS(paste0(dir, name, "/params.Rds"))
  total  = params$total
  m0     = params$m0
  A      = params$A
  rho_up = params$rho_up
  method = params$method
  impute_lb = params$impute_lb
  impute_ub = params$impute_ub
  start_batch = params$batch
  num_region = length(Y)
  num_island = length(beta)

  plots = output = vector("list", length(inits))
  names(plots) = names(output) = par_up = names(inits)
  for (batch in start_batch:60) {
    time = format(Sys.time(), "%a %b %d %X")
    cat("Batch", paste0(batch, ","), "Iteration", paste0(total, ","), time, "\r")
    T_inc = 100
    output$theta = array(dim = c(num_region, T_inc / 10))
    output$beta  = array(dim = c(num_island, T_inc / 10))
    output$sig2  = array(dim = c(T_inc / 10)            )
    output$tau2  = array(dim = c(T_inc / 10)            )
    output$Z     = array(dim = c(num_region, T_inc / 10))

    # Metropolis for Yikt
    t_accept = ifelse(t_accept < 1 / 6, 1 / 6, ifelse(t_accept > 0.75, 0.75, t_accept))
    theta_sd = ifelse(
      t_accept > 0.5,
      theta_sd * t_accept / 0.5,
      ifelse(t_accept < 0.35, theta_sd * t_accept / 0.35, theta_sd)
    )
    t_accept = rep(0, num_region)
    for(it in 1:T_inc) {
      #### impute missing Y's ####
      if (length(miss)) {
        if (method == "binom") {
          rate = expit(theta[miss])
          rp = stats::runif(
            length(miss),
            stats::pbinom(impute_lb - 0.1, round(n[miss]), rate),
            stats::pbinom(impute_ub + 0.1, round(n[miss]), rate)
          )
          Y[miss] = stats::qbinom(rp, round(n[miss]), rate)
        }
        if (method == "pois") {
          rate = n[miss] * exp(theta[miss])
          rp = stats::runif(
            length(miss),
            stats::ppois(impute_lb - 0.1, rate),
            stats::ppois(impute_ub + 0.1, rate)
          )
          Y[miss] = stats::qpois(rp, rate)
        }
      }

      Z = u_update_Z(Z, sig2, theta, beta, tau2, adjacency, num_adj, island_region, island_id)
      sig2 = u_update_sig2(sig2, Z, beta, tau2, adjacency, num_adj, island_region, num_island_region, A, m0, sig_a, sig_b)
      tau2 = u_update_tau2(tau2, theta, beta, Z, sig2, num_island_region, island_id, A, m0, tau_a, tau_b)
      theta = u_update_theta(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method)
      beta = u_update_beta(beta, theta, Z, tau2, sig2, A, m0, island_region)

      #### Save outputs ####
      if (it %% 10 == 0) {
        output$beta [, it / 10] = beta
        output$sig2 [it / 10]   = sig2
        output$tau2 [it / 10]   = tau2
        output$theta[, it / 10] = theta
        output$Z    [, it / 10] = Z
      }
      cat("Batch", paste0(batch, ","), "Iteration", paste0(total + it, ","), time, ifelse(it == T_inc, "\n", "\r"))
    }

    # modify meta-parameters, save outputs to respective files
    total = total + T_inc
    t_accept = t_accept / T_inc
    inits = list(
      theta = theta,
      beta = beta,
      Z = Z,
      sig2 = sig2,
      tau2 = tau2
    )
    priors$theta_sd = theta_sd
    priors$t_accept = t_accept
    params$total = total
    params$batch = batch
    saveRDS(params, paste0(dir, name, "/params.Rds"))
    saveRDS(priors, paste0(dir, name, "/priors.Rds"))
    saveRDS(inits,  paste0(dir, name, "/inits.Rds"))
    saveRDS(output$beta , paste0(dir, name, "/beta/" , "beta_out_" , batch, ".Rds"))
    saveRDS(output$theta, paste0(dir, name, "/theta/", "theta_out_", batch, ".Rds"))
    saveRDS(output$Z    , paste0(dir, name, "/Z/"    , "Z_out_"    , batch, ".Rds"))
    saveRDS(output$G    , paste0(dir, name, "/sig2/" , "sig2_out_" , batch, ".Rds"))
    saveRDS(output$tau2 , paste0(dir, name, "/tau2/" , "tau2_out_" , batch, ".Rds"))
    if (.show_plots) {
      # Output some of the estimates for plotting purposes
      plots$beta  = c(plots$beta,  output$beta [1, ])
      plots$theta = c(plots$theta, output$theta[1, ])
      plots$Z     = c(plots$Z,     output$Z    [1, ])
      plots$tau2  = c(plots$tau2,  output$tau2      )
      plots$sig2  = c(plots$sig2,  output$sig2      )

      grid = c(2, 3)
      graphics::par(mfrow = grid)
      burn = min(floor(total / 20), 200)
      its  = burn:(total / 10)
      plot(its * 10, plots$theta[its], type = "l", main = "theta")
      plot(its * 10, plots$beta[its], type = "l", main = "beta")
      plot(its * 10, plots$tau2[its], type = "l", main = "tau2")
      plot(its * 10, plots$sig2[its], type = "l", main = "sig2")
      plot(its * 10, plots$Z[its], type = "l", main = "Z")
    }

  }
  cat("Finished running model at:", format(Sys.time(), "%a %b %d %X"))
}
