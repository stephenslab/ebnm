# @title Sample from the posterior under point-normal prior
#
# @inheritParams ebnm_point_normal
#
# @param w The weight of the non-null component in the fitted prior.
#
# @param a The variance of the non-null component in the fitted prior.
#
# @param mu The mean of the distributions
#
# @param nsamp The number of samples to return per observation.
#
# @return An nsamp by length(x) matrix containing samples from the
#   posterior, with each row corresponding to a single sample.
#
#' @importFrom stats rbinom rnorm
#'
post_sampler_point_normal <- function(x, s, w, a, mu, nsamp) {
  wpost <- wpost_normal(x, s, w, a, mu)
  pmean_cond <- pmean_cond_normal(x, s, a, mu)
  pvar_cond <- pvar_cond_normal(s, a)

  is_nonnull <- rbinom(nsamp*length(wpost), 1, rep(wpost, nsamp))
  samp <- is_nonnull * rnorm(nsamp*length(wpost), rep(pmean_cond, nsamp),
                             rep(sqrt(pvar_cond), nsamp)) + (1 - is_nonnull) * mu
  return(matrix(samp, nrow=nsamp, byrow=T))
}
