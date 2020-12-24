#' @importFrom ashr my_etruncnorm my_e2truncnorm
#'
summary_results_point_laplace = function(x, s, w, a, mu, output) {
  x <- x - mu

  wpost <- wpost_laplace(x, s, w, a)
  lm <- lambda(x, s, a)

  post <- list()

  if (result_in_output(output)) {
    post$mean  <- wpost * (lm * my_etruncnorm(0, Inf, x - s^2 * a, s)
                           + (1 - lm) * my_etruncnorm(-Inf, 0, x + s^2 * a, s))
    post$mean2 <- wpost * (lm * my_e2truncnorm(0, Inf, x - s^2 * a, s)
                           + (1 - lm) * my_e2truncnorm(-Inf, 0, x + s^2 * a, s))
    if (any(is.infinite(s))) {
      post$mean[is.infinite(s)]  <- 0
      post$mean2[is.infinite(s)] <- 2 * w / a^2
    }
    post$sd <- sqrt(post$mean2 - post$mean^2)

    post$mean2 <- post$mean2 + mu^2 + 2 * mu * post$mean
    post$mean  <- post$mean + mu
  }

  if ("lfsr" %in% output) {
    post$lfsr <- (1 - wpost) + wpost * pmin(lm, 1 - lm)
    if (any(is.infinite(s))) {
      post$lfsr[is.infinite(s)] <- 1 - w / 2
    }
  }

  return(post)
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
#' @importFrom stats pnorm
#'
lambda <- function(x, s, a) {
  lm1 <- -a * x + pnorm(x / s - s * a, log.p = TRUE)
  lm2 <-  a * x + pnorm(x / s + s * a, log.p = TRUE, lower.tail = FALSE)

  lm <- 1 / (1 + exp(lm2 - lm1))

  return(lm)
}
