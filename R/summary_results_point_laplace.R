#' @importFrom ashr my_etruncnorm
#' @importFrom ashr my_e2truncnorm
#'
summary_results_point_laplace = function(x, s, w, a) {
  wpost <- wpost_laplace(x, s, w, a)
  lm <- lambda(x, s, a)

  posterior_mean <- wpost * (lm * my_etruncnorm(0, Inf, x - s^2 * a, s)
                             + (1 - lm) * my_etruncnorm(-Inf, 0, x + s^2 * a, s))
  if (any(is.infinite(s))) {
    posterior_mean[is.infinite(s)] <- 0
  }

  posterior_mean2 <- wpost * (lm * my_e2truncnorm(0, Inf, x - s^2 * a, s)
                              + (1 - lm) * my_e2truncnorm(-Inf, 0, x + s^2 * a, s))
  if (any(is.infinite(s))) {
    posterior_mean2[is.infinite(s)] <- 2 * w / a^2
  }

  return(data.frame(posterior_mean  = posterior_mean,
                    posterior_mean2 = posterior_mean2))
}

#  Calculate posterior weights for non-null effects.
#
#' @importFrom stats dnorm
#'
wpost_laplace <- function(x, s, w, a) {
  if (w == 0) {
    return(rep(0, length(x)))
  }

  if (w == 1) {
    return(rep(1, length(x)))
  }

  lf <- dnorm(x, 0, s, log = TRUE)
  lg <- logg_laplace(x, s, a)
  wpost <- w / (w + (1 - w) * exp(lf - lg))

  return(wpost)
}

# Compute the lambda function equation (2.7) from Kan Xu's thesis, which is
#   the posterior probability of being positive given a non-zero effect.
#
lambda <- function(x, s, a) {
  lm1 <- -a * x + pnorm(x / s - s * a, log.p = TRUE)
  lm2 <-  a * x + pnorm(x / s + s * a, log.p = TRUE, lower.tail = FALSE)

  lm <- 1 / (1 + exp(lm2 - lm1))

  return(lm)
}
