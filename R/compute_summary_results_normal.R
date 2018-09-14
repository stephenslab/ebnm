compute_summary_results_normal = function(x, s, w, a){
  wpost <- wpost_normal(x, s, w, a)
  pmean_cond <- pmean_cond_normal(x, s, a)
  pvar_cond <- pvar_cond_normal(s, a)

  PosteriorMean <- wpost * pmean_cond
  PosteriorMean2 <- wpost * (pmean_cond^2 + pvar_cond)

  return(data.frame(PosteriorMean = PosteriorMean,
                    PosteriorMean2 = PosteriorMean2))
}

#
#  Calculate the posterior weight for non-zero effect
#
#' @importFrom stats dnorm
#'
wpost_normal <- function(x, s, w, a) {
  if (w == 0) {
    return(rep(0, length(x)))
  }
  if (w == 1) {
    return(rep(1, length(x)))
  }

  lg <- dnorm(x, 0, sqrt(s^2 + 1/a), log = TRUE)
  lf <- dnorm(x, 0, s, log = TRUE)
  wpost <- w / (w + (1 - w) * exp(lf - lg))

  # Deal with zero sds:
  wpost[s == 0 & x == 0] <- 0
  wpost[s == 0 & x != 0] <- 1

  return(wpost)
}

#
#  Calculate the posterior mean for non-zero effect
#
pmean_cond_normal <- function(x, s, a) {
  return(x / (1 + s^2 * a))
}

#
#  Calculate the posterior variance for non-zero effect
#
pvar_cond_normal <- function(s, a) {
  return(s^2 / (1 + s^2 * a))
}
