# @title Sample from the posterior under normal prior
#
# @inheritParams ebnm_normal
#
# @param mu The mean of the normal in the fitted prior.
#
# @param a The precision of normal in the fitted prior.
#
# @param nsamp The number of samples to return per observation.
#
# @return An nsamp by length(x) matrix containing samples from the
#   posterior, with each row corresponding to a single sample.
#
#' @importFrom stats rbinom rnorm
#'
post_sampler_normal <- function(x, s, mu, a, nsamp) {
  pmean_cond <- pmean_cond_normal_mu(x, s, mu, a)
  pvar_cond <- pvar_cond_normal_mu(s, a)
  
  samp <- rnorm(nsamp*length(pmean_cond), rep(pmean_cond, nsamp),
                             rep(sqrt(pvar_cond), nsamp))
  return(matrix(samp, nrow = nsamp, byrow = T))
}
