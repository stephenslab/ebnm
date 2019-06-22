#' @describeIn ebnm Solves the EBNM problem using an ash prior.
#'
#' @importFrom ashr ash
#'
#' @export
#'
ebnm_ash = function(x,
                    s = 1,
                    mode = 0,
                    scale = "estimate",
                    g_init = NULL,
                    fix_g = FALSE,
                    output = output_default(),
                    ...) {
  ash_output <- output

  if ("result" %in% output) {
    ash_output <- setdiff(ash_output, "result")
    ash_output <- c(ash_output, "PosteriorMean", "PosteriorSD")
  }

  # Ash will accept either mode and mixsd or g, but not both.
  if (is.null(g_init)) {
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
    if ((!missing(mode) && !is.null(mode))
        || (!missing(scale) && !is.null(scale))) {
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
