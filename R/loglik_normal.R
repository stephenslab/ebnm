# Functions to compute the log likelihood under the normal prior with no point-mass and mean mu.
# a is the precision, like everywhere else

#' @importFrom stats dnorm
loglik_normal = function(x, s, mu, a) {
  sum(dnorm(x, mu, sqrt(s^2 + 1/a), log = TRUE))
}
