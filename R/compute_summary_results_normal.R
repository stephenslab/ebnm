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
  
  s2 = s * rep(1, length(x)) # new variable which properly recycles s with correct length
  
  reg_s_ind = is.finite(s2) & (s2 > 0)
  pm = x # initialize, case s_j=0
  pm[is.infinite(s2)] = mu # case s_j=Inf
  pm[reg_s_ind] = (x[reg_s_ind] + (s2[reg_s_ind]^2) * a * mu) / (1 + (s2[reg_s_ind]^2) * a)

  return(pm)
}

#
#  Calculate the posterior variance for non-zero effect
#  Note: don't need to worry about vector recycling here, taken care of when forming the dataframe above
#
pvar_cond_normal_mu <- function(s, a) {
  pvar_cond <- rep(1 / a, length(s))
  pvar_cond[is.finite(s)] <- s[is.finite(s)]^2 / (1 + s[is.finite(s)]^2 * a)
  
  return(pvar_cond)
}
