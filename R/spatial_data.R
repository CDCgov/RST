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

#' Check spatial data
#'
#' @noRd
check_spatial_data = function(adjacency, num_region) {
  cat("Checking spatial data...\n")
  # Check for errors
  errout = NULL
  errct  = 0
  # Length of adjacency not num_region
  if (length(adjacency) != num_region) {
    errct  = errct + 1
    errtxt = paste(errct, ": Adjacency different length than data. Ensure length(adjacency) == dim(Y)[1]")
    errout = c(errout, errtxt)
  }
  # Regions have no neighbors
  if (any(sapply(adjacency, length) == 0)) {
    errct  = errct + 1
    errtxt = paste(errct, ": Some regions have no neighbors. Ensure all regions have at least 1 neighbor. Check vignette('RST-adj') for more information")
    errout = c(errout, errtxt)
  }
  if (errct) {
    stop(paste(errct, "error(s) found in list 'data':\n", paste(errout, collapse = "\n ")))
  }
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
