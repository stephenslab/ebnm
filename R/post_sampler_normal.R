post_sampler_normal(x, s, w, a, nsamp) {
  wpost <- wpost_normal(x, s, w, a)
  pmean_cond <- pmean_cond_normal(x, s, a)
  pvar_cond <- pvar_cond_normal(s, a)

  is_nonnull <- rbinom(nsamp*length(wpost), 1, rep(wpost, nsamp))
  samp <- is_nonnull * rnorm(nsamp*length(wpost), rep(pmean_cond, nsamp),
                             rep(sqrt(pvar_cond), nsamp))
  return(matrix(samp, nrow=nsamp, byrow=T))
}
