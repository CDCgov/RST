#' Run Gibbs sampler
#' @param name Name of model and corresponding folder
#' @param dir  Directory where model lives
#' @param .show_plots Show or hide traceplots as 
#'
#' @export
run_sampler = function(name, dir, .show_plots = TRUE) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  model = readRDS(paste0(dir, name, "/params.Rds"))$model
  if (model == "ucar") {
    gibbs_u(name, dir, .show_plots)
  }
  if (model == "mcar") {
    gibbs_m(name, dir, .show_plots)
  }
  if (model == "ustcar") {
    gibbs_ust(name, dir, .show_plots)
  }
  if (model == "mstcar") {
    gibbs_mst(name, dir, .show_plots)
  }
}
