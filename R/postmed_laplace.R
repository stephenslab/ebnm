#
#  Find the posterior median for the Laplace prior for
#   given x (observations), s (sd), w and a.
#
postmed.laplace <- function(x, s, w, a) {

  if(w==0){return(rep(0,length(x)))}

  # Work with the absolute value of x, and for x > 25 use the approximation
  #  to dnorm(x-a)*beta.laplace(x, a)
	sx <- sign(x)
	x <- abs(x)
	xma <- x/s - s*a

	#The following code is on the non-log scale, and so can be numerically unstable.
	#We keep it here to help document where the log-based code comes from
	#zz <- 1/a * (1/s)*dnorm(xma) * (1/w + beta.laplace(x, s, a))
	#zz[xma > 25] <- 1/2
	#	mucor <- qnorm(pmin(zz, 1))

	logzz = log(1/a) + log(1/s) + dnorm(xma,log=TRUE) + log_inverse_w_plus_beta(w,x,s,a)
	mucor <- qnorm(pmin(logzz,0),log.p=TRUE)

	muhat <- sx * pmax(0, xma - mucor) * s
	return(muhat)
}

# this computes the log of the (1/w+beta) function used in computations of postmed
# it may not be the most efficient way. But it does it stably
log_inverse_w_plus_beta = function(w,x,s,a){
  lg = logg.laplace(x,s,a)
  lf = dnorm(x,0,s,log=TRUE)
  lfac = pmax(lg,lf)
  lfac + log(w*exp(lg-lfac)+(1-w)*exp(lf-lfac)) - log(w) -lf
}
