#' @describeIn ebnm Defines the default return values.
#'
#' @export
#'
output_default <- function() {
  return(c("result", "fitted_g", "loglik"))
}

#' @describeIn ebnm Lists all valid return values.
#'
#' @export
#'
output_all <- function() {
  return(c("result", "fitted_g", "loglik", "post_sampler", "lfsr"))
}
