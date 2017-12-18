# Functions to compute the log likelihood under the laplace prior
# More numerically stable than the original approach in negloglik_laplace.R

# This is the log of g, Laplace(a) convolved with normal, eqn (2.2) in Kan Xu's MS paper
logg_laplace = function(x,s,a){
  lg1 = -a*x + pnorm((x-s^2*a)/s,log=TRUE)
  lg2 = a*x + pnorm((x+s^2*a)/s,lower.tail = FALSE,log=TRUE)
  lfac = pmax(lg1,lg2) # avoid numeric issues by removing a log factor
  log(0.5) + log(a) + (a^2 * s^2)/2 + lfac + log(exp(lg1-lfac) + exp(lg2-lfac))
}

# return log((1-w)f + wg) as a vector
# deal with case w=1 and w=0 separately for stability
vloglik_laplace = function(x,s,w,a){
  if(w==0){return(dnorm(x/s,log=TRUE))}
  lg = logg_laplace(x,s,a)
  if(w==1){return(lg)}

  lf = dnorm(x/s,log=TRUE)
  lfac = pmax(lg,lf)
  return(lfac + log((1-w)*exp(lf-lfac) + w*exp(lg-lfac)))
}

loglik_laplace = function(x,s,w,a){
  sum(vloglik_laplace(x,s,w,a))
}
