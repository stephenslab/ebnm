# Functions to compute the log likelihood under the point-laplace prior.

loglik_point_laplace = function(x, s, w, a) {
  return(sum(vloglik_point_laplace(x, s, w, a)))
}

# Return log((1 - w)f + wg) as a vector (deal with cases w = 1 and w = 0
#   separately for stability).
#
vloglik_point_laplace = function(x, s, w, a) {
  if (w <= 0) {
    return(dnorm(x, sd = s, log = TRUE))
  }

  lg <- logg_laplace(x, s, a)
  if (w == 1) {
    return(lg)
  }

  lf <- dnorm(x, sd = s, log = TRUE)
  lfac <- pmax(lg, lf)
  return(lfac + log((1 - w) * exp(lf - lfac) + w * exp(lg - lfac)))
}


# This is the log of g, Laplace(a) convolved with normal, eqn (2.2) in Kan Xu's
#   MS paper.
#
logg_laplace = function(x, s, a) {
  lg1 <- -a * x + pnorm((x - s^2 * a) / s, log.p = TRUE)
  lg2 <-  a * x + pnorm((x + s^2 * a) / s, log.p = TRUE, lower.tail = FALSE)
  lfac <- pmax(lg1, lg2)
  return(log(a / 2) + s^2 * a^2 / 2 + lfac + log(exp(lg1 - lfac) + exp(lg2 - lfac)))
}
