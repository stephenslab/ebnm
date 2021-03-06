#' @importFrom ashr ash
#'
ebnm_ash_workhorse <- function(x,
                               s,
                               mode,
                               scale,
                               g_init,
                               fix_g,
                               output,
                               call,
                               ...) {
  if("mixsd" %in% names(list(...))) {
    stop("Use parameter 'scale' instead of 'mixsd'.")
  }
  if ("outputlevel" %in% names(list(...))) {
    stop("Use parameter 'output' instead of 'outputlevel'.")
  }

  if (identical(scale, "estimate")) {
    # Some ashr settings have implications for the grid:
    use_ashr_grid <- any(c("gridmult", "pointmass", "method") %in% names(list(...)))
    if(!identical(mode, "estimate") && !use_ashr_grid) {
      scale <- default_scale(x, s, mode)[-1]
    } else {
      # Let ashr do the grid estimation.
      scale <- NULL
    }
  }

  # Ash will accept either mode and mixsd or g, but not both.
  if (is.null(g_init)) {
    ash_res <- ash(betahat = as.vector(x),
                   sebetahat = as.vector(s),
                   mode = mode,
                   mixsd = scale,
                   fixg = fix_g,
                   outputlevel = ash_output(output),
                   ...)
  } else {
    if (!is.null(call$mode) || !is.null(call$scale)) {
      warning("mode and scale parameters are ignored when g_init is supplied.")
    }
    ash_res <- ash(betahat = as.vector(x),
                   sebetahat = as.vector(s),
                   g = g_init,
                   fixg = fix_g,
                   outputlevel = ash_output(output),
                   ...)
  }

  retlist <- list()

  if (posterior_in_output(output)) {
    posterior <- list()

    if (result_in_output(output)) {
      posterior$mean  <- ash_res$result$PosteriorMean
      posterior$sd    <- ash_res$result$PosteriorSD
      posterior$mean2 <- posterior$mean^2 + posterior$sd^2
    }

    if (lfsr_in_output(output)) {
      posterior$lfsr  <- ash_res$result$lfsr
    }

    retlist <- add_posterior_to_retlist(retlist, posterior, output)
  }

  if (g_in_output(output)) {
    retlist <- add_g_to_retlist(retlist, ash_res$fitted_g)
  }

  if (llik_in_output(output)) {
    retlist <- add_llik_to_retlist(retlist, ash_res$loglik)
  }

  if (sampler_in_output(output)) {
    retlist <- add_sampler_to_retlist(retlist, ash_res$post_sampler)
  }

  return(retlist)
}
