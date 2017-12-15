#' @title Solve the Empirical Bayes Normal Means problem with point-laplace prior
#' @details Given vectors of data x, and standard errors s, solve EBNM problem with "point-laplace" prior
#' (i.e. a mixture of point mass at 0 and Laplace distribution).
#' That is the prior is pi0 \delta_0 + (1-pi0)DExp(a) where Dexp is the double exponential
#' (Laplace) distribution, and (pi0,a) are esimated by marginal maximum likelihood. The model and code are based on
#' EbayesThresh by Johnstone and Silverman, but simplified by removing thresholding, and
#' improved numerical stability.
#' @param x a vector of observations
#' @param s a vector of standard deviations (or scalar if all equal)
#' @return a list with elements result, fitted_g, and loglik
#' @examples
#' mu = c(rep(0,1000), rexp(1000)) # means
#' s = rgamma(2000,1,1) #standard errors
#' x = mu + rnorm(2000,0,s) # observations
#' x.ebnm = ebnm_point_laplace(x,s)
#' ashr::get_pm(x.ebnm) # posterior mean
#'
#' @export
ebnm_point_laplace <- function (x,s=1) {
  #could consider makign more stable this way? But might have to be careful with log-likelihood
  #m_sdev <- mean(s)
  #s <- s/m_sdev
  #x <- x/m_sdev

	mle <- wandafromx.mle(x, s)
	w=mle$w
	a=mle$a

	loglik = loglik.laplace(x,s,w,a)
	result = compute_summary_results(x,s,w,a)

	#postmean <- postmean * m_sdev
	#postmean2 <- postmean2 * m_sdev^2

	retlist <- list(result=result, fitted_g = list(pi0  = 1-w, a = a), loglik = loglik)


	return(retlist)
}
