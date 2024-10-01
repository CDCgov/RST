#' Logit transformation
#'
#' @param x the value to be logit-transformed
logit = function(x) {
  log(x / (1 - x))
}

#' Expit transformation
#'
#' @param x the value to be expit-transformed
expit = function(x) {
  1 / (1 + exp(-x))
}
