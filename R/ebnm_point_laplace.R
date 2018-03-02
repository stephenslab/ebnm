#' @title Solve the Empirical Bayes Normal Means problem with point-laplace prior
#' @description Paragraph-length description goes here.
#'
#' @details Given vectors of data x, and standard errors s, solve EBNM problem with "point-laplace" prior
#' (i.e. a mixture of point mass at 0 and Laplace distribution).
#' That is the prior is \eqn{pi0 \delta_0 + (1-pi0)DExp(a)} where Dexp is the double exponential
#' (Laplace) distribution, and (pi0,a) are esimated by marginal maximum likelihood. The model and code are based on
#' EbayesThresh by Johnstone and Silverman, but simplified by removing thresholding, and
#' improved numerical stability.
#' @param x a vector of observations
#' @param s a vector of standard deviations (or scalar if all equal)
#' @param g The prior distribution (list with elements pi0,a). Usually this is unspecified (NULL) and
#' estimated from the data. However, it can be used in conjuction with fixg=TRUE
#' to specify the g to use (e.g. useful in simulations to do computations with the "true" g).
#' Or, if g is specified but fixg=FALSE, the g specifies the initial value of g used before optimization.
#' @param fixg If TRUE, don't estimate g but use the specified g.
#'
#' @return a list with elements result, fitted_g, and loglik
#' @examples
#' mu = c(rep(0,1000), rexp(1000)) # means
#' s = rgamma(2000,1,1) #standard errors
#' x = mu + rnorm(2000,0,s) # observations
#' x.ebnm = ebnm_point_laplace(x,s)
#' ashr::get_pm(x.ebnm) # posterior mean
#'
#' @export
ebnm_point_laplace <- function (x, s=1, g=NULL, fixg=FALSE, output=NULL) {
  output = set_output(output)

  #could consider making more stable this way? But might have to be careful with log-likelihood
  #m_sdev <- mean(s)
  #s <- s/m_sdev
  #x <- x/m_sdev

  if(is.null(g) & fixg){stop("must specify g if fixg=TRUE")}
  if(!is.null(g) & !fixg){stop("option to intialize from g not yet implemented")}

  # Estimate g from data
  if(!fixg){
    g <- mle_laplace(x, s)
  }

	w <- 1 - g$pi0
	a <- g$a

	retlist <- list()

	# Compute return values
	if ("summary_results" %in% output) {
	  result <- compute_summary_results_laplace(x, s, w, a)
	  #postmean <- postmean * m_sdev
	  #postmean2 <- postmean2 * m_sdev^2
	  retlist <- c(retlist, list(result=result))
	}
	if ("fitted_g" %in% output) {
	  retlist <- c(retlist, list(fitted_g=g))
	}
	if ("loglik" %in% output) {
	  loglik <- loglik_laplace(x, s, w, a)
	  retlist <- c(retlist, list(loglik=loglik))
	}
	if ("post_sampler" %in% output) {
	  stop("Posterior sampler not yet implemented for point-laplace prior.")
	}

	return(retlist)
}
