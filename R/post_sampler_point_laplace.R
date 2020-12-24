# Sample from the posterior under point-Laplace prior.
#
# @param nsamp The number of samples to return per observation.
#
# @return An nsamp by length(x) matrix containing samples from the
#   posterior, with each row corresponding to a single sample.
#
#' @importFrom stats rbinom
#'
#' @importFrom truncnorm rtruncnorm
#'
post_sampler_point_laplace <- function(x, s, w, a, mu, nsamp) {
  x <- x - mu

  wpost <- wpost_laplace(x, s, w, a)
  l <- lambda(x, s, a)

  nobs <- length(wpost)

  is_nonnull <- rbinom(nsamp * nobs, 1, rep(wpost, each = nsamp))
  is_positive <- rbinom(nsamp * nobs, 1, rep(l, each = nsamp))

  if (length(s) == 1) {
    s <- rep(s, nobs)
  }

  # rtruncnorm is not vectorized.
  negative_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, -Inf, 0, mean, sd)
  }, mean = x + s^2 * a, sd = s)
  positive_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, 0, Inf, mean, sd)
  }, mean = x - s^2 * a, sd = s)

  samp <- matrix(0, nrow = nsamp, ncol = length(wpost))
  samp[is_nonnull & is_positive] <- positive_samp[is_nonnull & is_positive]
  samp[is_nonnull & !is_positive] <- negative_samp[is_nonnull & !is_positive]

  samp <- samp + mu

  return(samp)
}
