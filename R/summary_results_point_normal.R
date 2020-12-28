pn_summres <- function(x, s, g, output) {
  pi0 <- g$pi0
  w   <- 1 - pi0
  a   <- g$a
  mu  <- g$mu

  wpost <- wpost_normal(x, s, w, a, mu)
  pmean_cond <- pmean_cond_normal(x, s, a, mu)
  pvar_cond <- pvar_cond_normal(s, a)

  posterior <- list()

  if (result_in_output(output)) {
    posterior$mean  <- wpost * pmean_cond + (1 - wpost) * mu
    posterior$mean2 <- wpost * (pmean_cond^2 + pvar_cond) + (1 - wpost) * (mu^2)
    posterior$sd    <- sqrt(posterior$mean2 - posterior$mean^2)
  }

  if (lfsr_in_output(output)) {
    posterior$lfsr  <- (1 - wpost) + wpost * pnorm(0, abs(pmean_cond), sqrt(pvar_cond))
  }

  return(posterior)
}

#  Calculate posterior weights for non-null effects.
#
wpost_normal <- function(x, s, w, a, mu) {
  if (w == 0) {
    return(rep(0, length(x)))
  }

  if (w == 1) {
    return(rep(1, length(x)))
  }

  llik.diff <- 0.5 * log(1 + 1 / (a * s^2))
  llik.diff <- llik.diff - 0.5 * (x - mu)^2 / (s^2 * (a * s^2 + 1))
  wpost <- w / (w + (1 - w) * exp(llik.diff))

  if (any(s == 0)) {
    wpost[s == 0 & x == mu] <- 0
    wpost[s == 0 & x != mu] <- 1
  }

  if (any(is.infinite(s))) {
    wpost[is.infinite(s)] <- w
  }

  return(wpost)
}

# Calculate posterior means for non-null effects.
#
pmean_cond_normal <- function(x, s, a, mu) {
  pm <- (x + s^2 * a * mu) / (1 + s^2 * a)

  if (any(is.infinite(s))) {
    pm[is.infinite(s)] <- mu
  }

  return(pm)
}

# Calculate posterior variances for non-null effects.
#
pvar_cond_normal <- function(s, a) {
  pvar_cond <- s^2 / (1 + s^2 * a)

  if (any(is.infinite(s))) {
    pvar_cond[is.infinite(s)] <- 1 / a
  }

  return(pvar_cond)
}

pn_partog <- function(par) {
  pi0 <- g$pi0
  w   <- 1 - pi0
  a   <- g$a
  mu  <- g$mu

  if (pi0 == 0) {
    g <- normalmix(pi = 1, mean = mu, sd = sqrt(1 / a))
  } else {
    g <- normalmix(pi = c(pi0, w),
                   mean = rep(mu, 2),
                   sd = c(0, sqrt(1 / a)))
  }

  return(g)
}
