# pi0 = 1, so G is the family of point masses \delta_0(\mu).
#
mle_point_only <- function(x, s, g, fix_a, fix_mu) {
  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }

  if (!fix_a) {
    # The value of a doesn't matter.
    g$a <- 1
  }

  if (!fix_mu) {
    if (length(s) == 1) {
      g$mu <- mean(x)
    } else {
      g$mu <- sum(x / s^2) / sum(1 / s^2)
    }
  }

  g$val <- loglik_point_normal(x, s, 0, g$a, g$mu)

  return(g)
}

# pi0 = 0, so G is the family of normal distributions N(\mu, 1 / a).
#
#' @importFrom stats optimize
#'
mle_normal <- function(x, s, g, control, fix_a, fix_mu) {
  if (fix_a) {
    # If a is fixed, the problem is equivalent to the "point only" problem.
    g <- mle_point_only(x, sqrt(s^2 + 1 / g$a), g, fix_a, fix_mu)
  } else if (length(s) == 1) {
    # If all SEs are identical, the solution has a simple closed form.
    if (!fix_mu) {
      g$mu <- mean(x)
    }
    g$a <- 1 / max(mean((x - g$mu)^2) - s^2, 0)
    g$val <- loglik_point_normal(x, s, 1, g$a, g$mu)
  } else {
    # Otherwise we need to solve a one-dimensional optimization problem.
    s2 <- s^2
    upper <- (max(x) - min(x))^2
    if (fix_mu) {
      # Only a needs to be estimated.
      xc2 <- (x - g$mu)^2
      optargs <- list(f = norm_llik, xc2 = xc2, s2 = s2, interval = c(0, upper),
                      maximum = TRUE)
      optres <- do.call(optimize, c(optargs, control))
    } else {
      # Both a and mu need to be estimated, but there is a closed-form
      #   expression for the optimal mu given a.
      optargs <- list(f = estmu_llik, x = x, s2 = s2, interval = c(0, upper),
                      maximum = TRUE)
      optres <- do.call(optimize, c(optargs, control))
      g$mu <- estmu(optres$maximum, x, s2)
    }
    g$a <- 1 / optres$maximum
    g$val <- optres$objective - length(x) * log(2 * pi) / 2
  }

  return(g)
}

# Used to estimate a when mu is fixed.
#
norm_llik <- function(sigma2, xc2, s2) {
  return(-0.5 * (sum(log(s2 + sigma2)) + sum(xc2 / (s2 + sigma2))))
}

# Used to estimate a when mu is not fixed.
#
estmu_llik <- function(sigma2, x, s2) {
  mu <- estmu(sigma2, x, s2)
  xc2 <- (x - mu)^2
  return(norm_llik(sigma2, xc2, s2))
}

# Given sigma2 = 1 / a, the optimal mu is:
#
estmu <- function(sigma2, x, s2) {
  return(sum(x / (s2 + sigma2)) / sum(1 / (s2 + sigma2)))
}
