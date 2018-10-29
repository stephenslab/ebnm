compute_summary_results_point_normal = function(x, s, w, a, mu){
  wpost <- wpost_normal(x, s, w, a, mu)
  pmean_cond <- pmean_cond_normal(x, s, a, mu)
  pvar_cond <- pvar_cond_normal(s, a)

  PosteriorMean <- wpost * pmean_cond + (1 - wpost) * mu
  PosteriorMean2 <- wpost * (pmean_cond^2 + pvar_cond) + (1 - wpost) * (mu^2)

  return(data.frame(PosteriorMean = PosteriorMean,
                    PosteriorMean2 = PosteriorMean2))
}

#
#  Calculate the posterior weight for non-null effect
#
#' @importFrom stats dnorm
#'
wpost_normal <- function(x, s, w, a, mu) {
  if (w == 0) {
    return(rep(0, length(x)))
  }
  if (w == 1) {
    return(rep(1, length(x)))
  }

  lg <- dnorm(x, mu, sqrt(s^2 + 1/a), log = TRUE)
  lf <- dnorm(x, mu, s, log = TRUE)
  wpost <- w / (w + (1 - w) * exp(lf - lg))

  # Deal with zero and infinite sds:
  wpost[s == 0 & x == mu] <- 0
  wpost[s == 0 & x != mu] <- 1
  wpost[is.infinite(s)] <- w

  return(wpost)
}

#
#  Calculate the posterior mean for non-null effect
#
pmean_cond_normal <- function(x, s, a, mu) {
  if (is.infinite(a)) { # if prior precision is infinite, posterior is prior mean (if s_j=0, still defer to prior)
    if (any(s == 0)) { # if prior precision infinite, but some s_j=0
      warning("Prior precision found to be infinite, but at least one s_j=0. Deferring to prior in this case")
    }
    return(rep(mu, length(x)))
  }
  
  pm = (x + (s^2) * a * mu) / (1 + (s^2) * a)
  pm[is.infinite(s)] = mu # case s_j=Inf
  
  return(pm)
}

#
#  Calculate the posterior variance for non-null effect
#
pvar_cond_normal <- function(s, a) {
  pvar_cond <- s^2 / (1 + s^2 * a)
  pvar_cond[is.infinite(s)] <- 1 / a

  return(pvar_cond)
}
