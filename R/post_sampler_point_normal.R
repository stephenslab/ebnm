#' # Sample from the posterior under point-normal prior
#' #
#' # @param nsamp The number of samples to return per observation.
#' #
#' # @return An nsamp by length(x) matrix containing samples from the
#' #   posterior, with each row corresponding to a single sample.
#' #
#' #' @importFrom stats rbinom rnorm
#' #'
#' post_sampler_point_normal <- function(x, s, w, a, mu, nsamp) {
#'   wpost <- wpost_normal(x, s, w, a, mu)
#'   pmean_cond <- pmean_cond_normal(x, s, a, mu)
#'   pvar_cond <- pvar_cond_normal(s, a)
#'
#'   nobs <- length(x)
#'   is_nonnull <- rbinom(nsamp * nobs, 1, rep(wpost, each = nsamp))
#'   samp <- is_nonnull * rnorm(nsamp * nobs,
#'                              mean = rep(pmean_cond, each = nsamp),
#'                              sd = rep(sqrt(pvar_cond), each = nsamp))
#'   samp <- samp + (1 - is_nonnull) * mu
#'
#'   return(matrix(samp, nrow = nsamp))
#' }
