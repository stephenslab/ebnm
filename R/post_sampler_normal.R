#' @title Sample from the posterior for data (x,s) under point normal prior
#' @param x observations
#' @param s standard deviations
#' @param w weight of non-null component in fitted prior
#' @param a variance of non-null component in fitted prior
#' @param nsamp number of samples to return per observation
post_sampler_normal(x, s, w, a, nsamp) {
  wpost <- wpost_normal(x, s, w, a)
  pmean_cond <- pmean_cond_normal(x, s, a)
  pvar_cond <- pvar_cond_normal(s, a)

  is_nonnull <- rbinom(nsamp*length(wpost), 1, rep(wpost, nsamp))
  samp <- is_nonnull * rnorm(nsamp*length(wpost), rep(pmean_cond, nsamp),
                             rep(sqrt(pvar_cond), nsamp))
  return(matrix(samp, nrow=nsamp, byrow=T))
}
