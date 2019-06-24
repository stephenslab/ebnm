#' @describeIn ebnm A generic function for solving the EBNM problem using an
#'   adaptive shrinkage (ash) prior.
#'
#' @export
#'
ebnm_ash <- function(x,
                     s = 1,
                     mode = 0,
                     scale = "estimate",
                     g_init = NULL,
                     fix_g = FALSE,
                     output = output_default(),
                     ...) {
  return(ebnm_ash_workhorse(x = x,
                            s = s,
                            mode = mode,
                            scale = scale,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            call = match.call(),
                            ...))
}

# The workhorse function is used by all ebnm functions that call into ashr.
#
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
  ash_output <- output

  if ("result" %in% output) {
    ash_output <- setdiff(ash_output, "result")
    ash_output <- c(ash_output, "PosteriorMean", "PosteriorSD")
  }

  # Ash will accept either mode and mixsd or g, but not both.
  if (is.null(g_init)) {
    # Allow partial matching for mode and scale.
    if (identical(pmatch(mode, "estimate"), 1L)) {
      mode <- "estimate"
    }
    if (identical(pmatch(scale, "estimate"), 1L)) {
      scale <- NULL
    }
    ash.res <- ash(betahat = as.vector(x),
                   sebetahat = as.vector(s),
                   mode = mode,
                   mixsd = scale,
                   fixg = fix_g,
                   outputlevel = ash_output,
                   ...)
  } else {
    if (!is.null(call$mode) || !is.null(call$scale)) {
      warning("mode and scale parameters are ignored when g_init is supplied.")
    }
    ash.res <- ash(betahat = as.vector(x),
                   sebetahat = as.vector(s),
                   g = g_init,
                   fixg = fix_g,
                   outputlevel = ash_output,
                   ...)
  }

  res <- list()

  if ("result" %in% output || "lfsr" %in% output) {
    res$result <- list()
    if ("result" %in% output) {
      pm <- ash.res$result$PosteriorMean
      res$result$posterior_mean <- pm
      res$result$posterior_mean2 <- pm^2 + ash.res$result$PosteriorSD^2
    }
    if ("lfsr" %in% output) {
      res$result$lfsr <- ash.res$result$lfsr
    }
    res$result <- data.frame(res$result)
  }

  if ("fitted_g" %in% output) {
    res$fitted_g <- ash.res$fitted_g
  }

  if ("loglik" %in% output) {
    res$loglik <- ash.res$loglik
  }

  if ("post_sampler" %in% output) {
    res$post_sampler <- ash.res$post_sampler
  }

  return(res)
}
