#' @importFrom ashr normalmix
#'
flat_workhorse <- function(x,
                           s,
                           mode,
                           scale,
                           g_init,
                           fix_g,
                           output,
                           call,
                           ...) {
  retlist <- list()

  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, x, s)
  }

  if (posterior_in_output(output)) {
    posterior <- list()

    if (result_in_output(output)) {
      posterior$mean  <- x
      posterior$sd    <- s
      posterior$mean2 <- x^2 + s^2
    }

    if (lfsr_in_output(output)) {
      posterior$lfsr  <- pnorm(-abs(x) / s)
    }

    retlist <- add_posterior_to_retlist(retlist, posterior, output, x)
  }

  if (g_in_output(output)) {
    fitted_g <- ashr::normalmix(1, 0, Inf)
    retlist <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    # Log likelihood is not well defined so just return zero.
    retlist <- add_llik_to_retlist(retlist, NA, x, df = 0)
  }

  if (sampler_in_output(output)) {
    post_sampler <- function(nsamp) {
      samp <- matrix(
        rnorm(nsamp * length(x), x, s),
        nrow = nsamp,
        byrow = TRUE
      )
      colnames(samp) <- names(x)
      return(samp)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(retlist)
}
