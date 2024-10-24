#' Gibbs sampler
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#'
#' @noRd
gibbs_mst = function(name, dir, .show_plots) {
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

  priors = readRDS(paste0(dir, name, "/priors.Rds"))
  Ag_scale = priors$Ag_scale
  G_df     = priors$G_df
  Ag_df    = priors$Ag_df
  tau_a    = priors$tau_a
  tau_b    = priors$tau_b
  rho_a    = priors$rho_a
  rho_b    = priors$rho_b
  theta_sd = priors$theta_sd
  rho_sd   = priors$rho_sd
  t_accept = priors$t_accept
  r_accept = priors$r_accept

  inits = readRDS(paste0(dir, name, "/inits.Rds"))
  theta = inits$theta
  beta  = inits$beta
  Z     = inits$Z
  G     = inits$G
  rho   = inits$rho
  tau2  = inits$tau2
  Ag    = inits$Ag

  params = readRDS(paste0(dir, name, "/params.Rds"))
  total  = params$total
  rho_up = params$rho_up
  method = params$method
  impute_lb = params$impute_lb
  impute_ub = params$impute_ub
  start_batch = params$batch

  plots = output = vector("list", length(inits))
  names(plots) = names(output) = par_up = names(inits)
  if (!rho_up) par_up = par_up[-which(par_up == "rho")]
  for (batch in start_batch:60) {
    time = format(Sys.time(), "%a %b %d %X")
    cat("Batch", paste0(batch, ","), "Iteration", paste0(total, ","), time, "\r")
    T_inc = 100
    output$theta = array(dim = c(dim(theta)  , T_inc / 10))
    output$beta  = array(dim = c(dim(beta)   , T_inc / 10))
    output$G     = array(dim = c(dim(G)      , T_inc / 10))
    output$tau2  = array(dim = c(length(tau2), T_inc / 10))
    output$Ag    = array(dim = c(dim(Ag)     , T_inc / 10))
    output$Z     = array(dim = c(dim(Z)      , T_inc / 10))
    if (rho_up) {
      output$rho = array(dim = c(length(rho) , T_inc / 10))
    }

    r_accept = ifelse(r_accept < 1 / 6, 1 / 6, ifelse(r_accept > 0.75, 0.75, r_accept))
    rho_sd = ifelse(
      r_accept > 0.5,
      rho_sd * r_accept / 0.5,
      ifelse(r_accept < 0.25, rho_sd * r_accept / 0.25, rho_sd)
    )
    r_accept = rep(0, length(rho))

    # Metropolis for Yikt
    t_accept = ifelse(t_accept < 1 / 6, 1 / 6, ifelse(t_accept > 0.75, 0.75, t_accept))
    theta_sd = ifelse(
      t_accept > 0.5,
      theta_sd * t_accept / 0.5,
      ifelse(t_accept < 0.35, theta_sd * t_accept / 0.35, theta_sd)
    )
    t_accept = array(0, dim = dim(theta))
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

      ##### update beta ####
      beta = mst_update_beta(beta, theta, Z, tau2, island_region)

      #### update Z ####
      Z = mst_update_Z(Z, G, theta, beta, rho, tau2, adjacency, num_adj, island_region, island_id)

      #### update G ####
      G = mst_update_G(G, Z, Ag, rho, G_df, adjacency, num_island)

      #### update Ag ####
      Ag = mst_update_Ag(Ag, G, Ag_scale, G_df, Ag_df)

      ##### update tau2 ####
      tau2 = mst_update_tau2(tau2, theta, beta, Z, tau_a, tau_b, island_id)

      #### update theta ####
      theta = mst_update_theta(theta, t_accept, Y, n, Z, beta, tau2, theta_sd, island_id, method)

      #### update rho ####
      if (rho_up) {
        rho = mst_update_rho(rho, r_accept, G, Z, rho_a, rho_b, rho_sd, adjacency, num_island)
      }

      #### Save outputs ####
      if (it %% 10 == 0) {
        output$beta [, , , it / 10] = beta
        output$G    [, , , it / 10] = G
        output$tau2 [,     it / 10] = tau2
        output$theta[, , , it / 10] = theta
        output$Z    [, , , it / 10] = Z
        output$Ag   [, ,   it / 10] = Ag
        if (rho_up) {
          output$rho[,     it / 10] = rho
        }
      }
      cat("Batch", paste0(batch, ","), "Iteration", paste0(total + it, ","), time, ifelse(it == T_inc, "\n", "\r"))
    }

    # modify meta-parameters, save outputs to respective files
    total = total + T_inc
    r_accept = r_accept / T_inc
    t_accept = t_accept / T_inc
    inits = list(
      theta = theta,
      beta = beta,
      Z = Z,
      G = G,
      rho = rho,
      tau2 = tau2,
      Ag = Ag
    )
    priors$theta_sd = theta_sd
    priors$rho_sd = rho_sd
    priors$t_accept = t_accept
    priors$r_accept = r_accept
    params$total = total
    params$batch = batch
    saveRDS(params, paste0(dir, name, "/params.Rds"))
    saveRDS(priors, paste0(dir, name, "/priors.Rds"))
    saveRDS(inits,  paste0(dir, name, "/inits.Rds"))
    saveRDS(output$beta , paste0(dir, name, "/beta/" , "beta_out_" , batch, ".Rds"))
    saveRDS(output$theta, paste0(dir, name, "/theta/", "theta_out_", batch, ".Rds"))
    saveRDS(output$Z    , paste0(dir, name, "/Z/"    , "Z_out_"    , batch, ".Rds"))
    saveRDS(output$G    , paste0(dir, name, "/G/"    , "G_out_"    , batch, ".Rds"))
    saveRDS(output$tau2 , paste0(dir, name, "/tau2/" , "tau2_out_" , batch, ".Rds"))
    saveRDS(output$Ag   , paste0(dir, name, "/Ag/"   , "Ag_out_"   , batch, ".Rds"))
    if (rho_up) {
      saveRDS(output$rho, paste0(dir, name, "/rho/", "rho_out_", batch, ".Rds"))
    }
    if (.show_plots) {
      # Output some of the estimates for plotting purposes
      plots$beta  = c(plots$beta,  output$beta [1, 1, 1, ])
      plots$theta = c(plots$theta, output$theta[1, 1, 1, ])
      plots$Z     = c(plots$Z,     output$Z    [1, 1, 1, ])
      plots$tau2  = c(plots$tau2,  output$tau2 [1,       ])
      plots$G     = c(plots$G,     output$G    [1, 1, 1, ])
      plots$Ag    = c(plots$Ag,    output$Ag   [1, 1,    ])
      if (rho_up) {
        plots$rho = c(plots$rho, output$rho[1, ])
      }

      grid = c(2, 3)
      if (rho_up) grid = c(2, 4)
      graphics::par(mfrow = grid)
      burn = min(floor(total / 20), 200)
      its  = burn:(total / 10)
      plot(its * 10, plots$theta[its], type = "l", main = "theta")
      plot(its * 10, plots$beta[its], type = "l", main = "beta")
      plot(its * 10, plots$tau2[its], type = "l", main = "tau2")
      plot(its * 10, plots$G[its], type = "l", main = "G")
      plot(its * 10, plots$Z[its], type = "l", main = "Z")
      plot(its * 10, plots$Ag[its], type = "l", main = "Ag")
      if (rho_up) {
        plot(its * 10, plots$rho[its], type = "l", main = "rho")
      }
    }

  }
  cat("Finished running model at:", format(Sys.time(), "%a %b %d %X"))
}
