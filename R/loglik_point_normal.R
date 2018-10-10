# Functions to compute the log likelihood under the normal prior.

loglik_point_normal = function(x, s, w, a) {
  sum(vloglik_point_normal(x, s, w, a))
}

# Return log((1 - w)f + wg) as a vector (deal with cases w = 1 and w = 0
#   separately for stability).
#
#' @importFrom stats dnorm
#'
vloglik_point_normal = function(x, s, w, a) {
  if (w <= 0) {
    return(dnorm(x, 0, s, log = TRUE))
  }

  lg <- dnorm(x, 0, sqrt(s^2 + 1/a), log = TRUE)
  if (w >= 1) {
    return(lg)
  }

  lf <- dnorm(x, 0, s, log = TRUE)
  lfac <- pmax(lg, lf)
  result <- lfac + log((1 - w) * exp(lf - lfac) + w * exp(lg - lfac))

  # Deal with zero sds:
  result[s == 0 & x == 0] <- log(1 - w)
  result[s == 0 & x != 0] <- log(w) + lg[s == 0 & x != 0]

  return(result)
}
