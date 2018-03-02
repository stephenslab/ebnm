set_output(output) {
  if (is.null(output)) {
    output <- c("summary_results", "fitted_g", "loglik")
  }
  return(output)
}
