#' Run Gibbs sampler
#' @param name Name of model and corresponding folder
#' @param dir  Directory where model lives
#'
#' @export
run_sampler = function(name, dir) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir = paste0(dir, "/")
  }
  gibbs_mst(name, dir)
}
