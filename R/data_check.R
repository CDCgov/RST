#' Check data
#'
#' @noRd
check_data = function(data) {
  cat("Checking data...\n")
  Y = data$Y
  n = data$n
  chk = c("Y", "n")
  miss = sapply(1:length(chk), \(x) !any(names(data) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'data': ", paste(chk[miss], collapse = ", "))
  }

  # Check for warnings
  warnout = NULL
  warnct  = 0
  # Check for unused elements in 'data'
  chk_elem = which(!(names(data) %in% c("Y", "n")))
  if (length(chk_elem)) {
    warnct  = warnct + 1
    warntxt = paste(warnct, ": Unused elements of list 'data':", paste(names(data)[chk_elem], collapse = ", "))
    warnout = c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'data':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout = NULL
  errct  = 0
  # Dimensions of Y and n are not the same
  dimtest = NULL
  if (is.null(dim(Y))) {
    dimtest = length(Y) != length(n)
  } else {
    dimtest = any(dim(Y) != dim(n))
  }
  if (dimtest) {
    errct  = errct + 1
    errtxt = paste(errct, ": Data not same dimensions. Ensure dim(Y) == dim(n)")
    errout = c(errout, errtxt)
  }
  # Values of Y are either negative or infinite
  Ychk = Y[which(!is.na(Y) & !is.null(Y))]
  if (any((Ychk < 0) | is.infinite(Ychk))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Invalid Y values. Check that all Y's are at least 0 and finite")
    errout = c(errout, errtxt)
  }
  # Values of n are either negative or infinite
  if (any((n < 0) | is.infinite(n))) {
    errct  = errct + 1
    errtxt = paste(errct, ": Invalid n values. Check that all n's are at least 0 and finite")
    errout = c(errout, errtxt)
  }
  if (errct) {
    stop(paste(errct, "error(s) found in list 'data':\n", paste(errout, collapse = "\n ")))
  }
}
