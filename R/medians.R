#' Calculate Medians
#' @param output Array of output to generate estimates for.
#'
#' @export
get_medians = function(output) {
  ndims = length(dim(output)) - 1
  apply(output, 1:ndims, stats::median)
}

#' Calculate Relative Precisions
#' @param output Array of output to generate relative precisions for.
#' @param medians Array of medians generated from output
#' @param perc_ci Number from 0 to 1. Generates relative precisions for a `perc_ci * 100`\% credible interval.
#'
#' @export
get_relative_precision = function(output, medians, perc_ci) {
  ndims = length(dim(output)) - 1
  alpha = (1 - perc_ci) / 2
  ci = apply(output, 1:ndims, stats::quantile, c(alpha, 1 - alpha))
  medians / (ci[2, , , ] - ci[1, , , ])
}
