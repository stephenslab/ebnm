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
  if (scale == "estimate") {
    scale <- NULL
  }

  ash_output <- output
  if ("result" %in% output) {
    ash_output <- setdiff(ash_output, "result")
    ash_output <- c(ash_output, "PosteriorMean", "PosteriorSD")
  }

  # Ash will accept either mode and mixsd or g, but not both.
  if (is.null(g_init)) {
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

  retlist <- list()

  if ("result" %in% output || "lfsr" %in% output) {
    retlist$result <- list()
    if ("result" %in% output) {
      pm  <- ash.res$result$PosteriorMean
      pm2 <- ash.res$result$PosteriorSD^2
      retlist$result$posterior_mean  <- pm
      retlist$result$posterior_mean2 <- pm^2 + pm2
    }
    if ("lfsr" %in% output) {
      retlist$result$lfsr <- ash.res$result$lfsr
    }
    retlist$result <- data.frame(retlist$result)
  }

  if ("fitted_g" %in% output) {
    retlist$fitted_g <- ash.res$fitted_g
  }

  if ("loglik" %in% output) {
    retlist$loglik <- ash.res$loglik
  }

  if ("post_sampler" %in% output) {
    retlist$post_sampler <- ash.res$post_sampler
  }

  return(retlist)
}
