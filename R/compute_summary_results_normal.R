compute_summary_results_normal = function(x, s, mu, a){
  pmean_cond <- pmean_cond_normal_mu(x, s, mu, a)
  pvar_cond <- pvar_cond_normal_mu(s, a)
  
  PosteriorMean <- pmean_cond
  PosteriorMean2 <- (pmean_cond^2 + pvar_cond)
  
  return(data.frame(PosteriorMean = PosteriorMean,
                    PosteriorMean2 = PosteriorMean2))
}



#
#  Calculate the posterior mean effect
#
pmean_cond_normal_mu <- function(x, s, mu, a) {
  if (is.infinite(a)) { # if prior precision is infinite, posterior is proir mean (if s_j=0, still defer to prior)
    return(rep(mu, length(x)))
  }
  
  return(sapply(1:length(x), function(l) ifelse(s[l] == 0, x[l], weighted.mean(c(x[l], mu), c(1/s[l]^2, a)))))
}

#
#  Calculate the posterior variance for non-zero effect
#
pvar_cond_normal_mu <- function(s, a) {
  pvar_cond <- rep(1 / a, length(s))
  pvar_cond[is.finite(s)] <- s[is.finite(s)]^2 / (1 + s[is.finite(s)]^2 * a)
  
  return(pvar_cond)
}
