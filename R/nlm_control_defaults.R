# Default control parameter for nlm (used by point_laplace and point_normal).

nlm_control_defaults <- function() {
  return(list(ndigit = 8, stepmax = 10, check.analyticals = FALSE))
}
