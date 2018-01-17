#' @title Solve the Empirical Bayes Normal Means problem with point-normal prior
#' @description Paragraph-length description goes here.
#'
#' @details Given vectors of data x, and standard errors s, solve EBNM problem with "point-normal" prior
#' (i.e. a mixture of point mass at 0 and normal distribution).
#' That is the prior is \eqn{pi0 \delta_0 + (1-pi0)N(0,1/a)} where N is the normal
#' distribution, and (pi0,a) are estimated by marginal maximum likelihood.
#' @param x a vector of observations
#' @param s a vector of standard deviations (or scalar if all equal)
#' @param g The prior distribution (list with elements pi0,a). Usually this is unspecified (NULL) and
#' estimated from the data. However, it can be used in conjuction with fixg=TRUE
#' to specify the g to use (e.g. useful in simulations to do computations with the "true" g).
#' Or, if g is specified but fixg=FALSE, the g specifies the initial value of g used before optimization.
#' @param fixg If TRUE, don't estimate g but use the specified g.
#' @param norm normalization factor to divide x and s by before running optimization (should not affect results, but improves numerical stability when x and s are tiny)
#'
#' @return a list with elements result, fitted_g, and loglik
#' @examples
#' mu = c(rep(0,1000), rexp(1000)) # means
#' s = rgamma(2000,1,1) #standard errors
#' x = mu + rnorm(2000,0,s) # observations
#' x.ebnm = ebnm_point_normal(x,s)
#' ashr::get_pm(x.ebnm) # posterior mean
#'
#' @export
ebnm_point_normal <- function (x,s=1,g=NULL,fixg=FALSE, norm = mean(s)) {
  #could consider makign more stable this way? But might have to be careful with log-likelihood

  s <- s/norm
  x <- x/norm

  if(is.null(g) & fixg){stop("must specify g if fixg=TRUE")}
  if(!is.null(g) & !fixg){stop("option to intialize from g not yet implemented")}

  # Estimate g from data
  if(!fixg){
    g <- mle_normal_logscale_grad(x, s)
  }

	w=1-g$pi0
	a=g$a

	loglik = loglik_normal(x,s,w,a)
	result = compute_summary_results_normal(x,s,w,a)

	# adjust results back to original scale
	result$PosteriorMean <- result$PosteriorMean * norm
	result$PosteriorMean2 <- result$PosteriorMean2 * norm^2
	g$a <- g$a / (norm^2)
	loglik = loglik - length(x)*log(norm)


	retlist <- list(result=result, fitted_g = g, loglik = loglik)


	return(retlist)
}
