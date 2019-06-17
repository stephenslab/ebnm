#' @describeIn ebnm Solve the EBNM problem using an ash prior.
#'
#' @export
#'
ebnm_ash = function(x,
                    s = 1,
                    g_init = NULL,
                    fix_g = FALSE,
                    output = output_default(),
                    ...) {
  ash_output <- output
  if ("result" %in% output) {
    ash_output <- setdiff(ash_output, "result")
    ash_output <- c(ash_output, "PosteriorMean", "PosteriorSD")
  }

  res <- ashr::ash(betahat = as.vector(x),
                   sebetahat = as.vector(s),
                   g = g_init,
                   fixg = fix_g,
                   outputlevel = ash_output,
                   ...)

  if ("result" %in% output) {
    res$result$betahat <- NULL
    res$result$sebetahat <- NULL
    psd <- res$result$PosteriorSD
    res$result$PosteriorSD <- NULL
    res$result$PosteriorMean2 <- res$result$PosteriorMean^2 + psd^2
  }

  if (!is.null(res$result)) {
    res$result <- as.list(res$result)
  }

  return(res)
}
