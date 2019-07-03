# Default arguments for nlm parameters (used by point_laplace, point_normal,
#   and normal).

nlm_control_defaults <- function() {
  return(list(ndigit = 8, stepmax = 10, check.analyticals = FALSE))
}
