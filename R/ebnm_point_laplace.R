#' @describeIn ebnm Solve the EBNM problem using a point-laplace prior.
#'
#' @export
#'
ebnm_point_laplace <- function (x, s=1, g=list(), fixg=FALSE, output=NULL) {
  output <- set_output(output)
  check_args(x, s, g, fixg, output)

  if (!fixg && (length(g) > 0)) {
    stop("Option to intialize from g not yet implemented for 'point_laplace' ",
         "priors.")
  }
  if ("post_sampler" %in% output) {
    stop("Posterior sampler not yet implemented for 'point-laplace' priors.")
  }

  #could consider making more stable this way? But might have to be careful with log-likelihood
  #m_sdev <- mean(s)
  #s <- s/m_sdev
  #x <- x/m_sdev

  # Estimate g from data
  if (!fixg) {
    g <- mle_laplace(x, s)
  }
  g$mu <- 0

	w <- 1 - g$pi0
	a <- g$a

	retlist <- list()

	# Compute return values
	if ("result" %in% output) {
	  result <- compute_summary_results_laplace(x, s, w, a)
	  #postmean <- postmean * m_sdev
	  #postmean2 <- postmean2 * m_sdev^2
	  retlist <- c(retlist, list(result = result))
	}
	if ("fitted_g" %in% output) {
	  retlist <- c(retlist, list(fitted_g = g))
	}
	if ("loglik" %in% output) {
	  loglik <- loglik_laplace(x, s, w, a)
	  retlist <- c(retlist, list(loglik = loglik))
	}

	return(retlist)
}
