#' Get spatial data
#'
#' @noRd
get_spatial_data = function(adjacency, num_region, .ignore_checks) {
  if (!.ignore_checks) {
    check_spatial_data(adjacency, num_region)
  }
  num_adj = sapply(adjacency, length)
  island_region = lapply(get_islands(adjacency), \(x) x - 1)
  adjacency = lapply(adjacency, \(x) x - 1)
  num_island = length(island_region)
  island_id = rep(NA, num_region)
  for (isl in 1:num_island) {
    island_id[island_region[[isl]] + 1] = isl - 1
  }
  list(
    adjacency = adjacency,
    num_adj = num_adj,
    island_region = island_region,
    island_id = island_id,
    num_island = num_island
  )
}

#' Get islands
#'
#' @noRd
get_islands = function(adjacency) {
  f = 1:length(adjacency)
  island_region = list(); group = 0
  while(length(f) > 0) {
    active_list   = f[1]
    inactive_list = NULL
    t = 0
    while(length(active_list) > 0) {
      Na = adjacency[[active_list[1]]]
      active_list   = unique(c(active_list, Na[which(!(Na %in% inactive_list))]))
      inactive_list = c(inactive_list, active_list[1])
      active_list   = active_list[-1]
    }
    group = group + 1
    inactive_list = sort(inactive_list)
    island_region[[group]] = inactive_list
    f = setdiff(f, inactive_list)
  }
  island_region
}
