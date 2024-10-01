#' Age-standardize samples
#' @param sample Sample array to age-standardize
#' @param std_pop A vector of standard populations
#' @param margin The margin on which the age data is stratified
#'
#' @export
age_standardize = function(sample, std_pop, margin) {
  mar = seq_along(dim(sample))[-margin]
  wts = std_pop / sum(std_pop)
  new_dim = dim(sample)
  new_dim[margin] = 1
  std_samp = array(
    apply(sweep(sample, margin, wts, "*"), mar, sum, na.rm = TRUE),
    dim = new_dim
  )
  std_samp
}

#' Aggregate samples
#' @param sample Sample array to age-standardize
#' @param pop A population array to be used for weighted averages
#' @param margin The margin on which the groups of interest are stratified
#'
#' @export
group_aggregate = function(sample, pop, margin) {
  mar = seq_along(dim(sample))[-margin]
  pop_arr = array(pop, dim = c(dim(pop), rev(dim(sample))[1]))
  new_dim = dim(sample)
  new_dim[margin] = 1
  agg_samp = array(
    apply(sample * pop_arr, mar, sum, na.rm = TRUE) / apply(pop_arr, mar, sum, na.rm = TRUE),
    dim = new_dim
  )
  agg_samp
}

#' Bind and name standardized/aggregated samples
#' @param sample Original sample array
#' @param agg_sample Newly aggregated/standardized array
#' @param margin The margin on which the groups of interest are stratified
#' @param new_name The name of the newly aggregated group
#'
#' @export
bind_samples = function(sample, agg_sample, margin, new_name = NULL) {
  if (!is.null(new_name)) {
    newnames = c(dimnames(sample)[[margin]], new_name)
  }
  sample = abind::abind(sample, agg_sample, along = margin)
  dimnames(sample)[[margin]] = newnames
  sample
}

#' Load population array
#' @param name Name of the model
#' @param dir Directory of the model
#'
#' @export
load_pop = function(name, dir) {
  readRDS(paste0(dir, "/", name, "/data.Rds"))$n
}

#' Aggregate population arrays
#' @param pop Population array
#' @param margin The margin on which the groups of interest are stratified
#'
#' @export
pop_aggregate = function(pop, margin) {
  mar = seq_along(dim(pop))[-margin]
  new_dim = dim(pop)
  new_dim[margin] = 1
  agg_samp = array(apply(pop, mar, sum, na.rm = TRUE), dim = new_dim)
  agg_samp
}
