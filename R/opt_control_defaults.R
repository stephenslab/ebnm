# Default arguments for optimization parameters (used by normal, point_normal,
#   point_laplace, and point_exponential).

nlm_control_defaults <- function() {
  # Sometimes nlm thinks that the gradient is being calculated incorrectly.
  #   Reducing the number of significant digits often solves the problem.
  return(list(ndigit = 8, stepmax = 5, check.analyticals = FALSE))
}

lbfgsb_control_defaults <- function() {
  return(list())
}

trust_control_defaults <- function () {
  return(list(rinit = 2, rmax = 5))
}

optimize_control_defaults <- function() {
  return(list())
}
