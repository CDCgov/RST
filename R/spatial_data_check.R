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
