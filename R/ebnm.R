#' @title Solve the Empirical Bayes Normal Means problem
#' @details Given vectors of data x, and standard errors s, fit the point laplace prior
#' (i.e. a mixture of point mass at 0 and Laplace distribution)
#' by fitting the two parameters (pi0 and laplace rate). The model and code are based on
#' EbayesThresh by Johnstone and Silverman, but simplified by removing thresholding, and
#' improved numerical stability
#' @param x a vector of observations
#' @param s a vector of standard deviations (or scalar if all equal)
#' @return a list with elements result, fitted_g, and loglik
#' @export
ebnm <- function (x,s=1) {
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
