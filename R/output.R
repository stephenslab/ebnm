#
# Sets output for ebnm_point_laplace and ebnm_point_normal. Defaults to
#   summary results, fitted prior, and log likelihood.
#
set_output <- function(output) {
  if (is.null(output)) {
    output <- c("summary_results", "fitted_g", "loglik")
  }
  return(output)
}
