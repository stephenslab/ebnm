#' @importFrom ashr my_etruncnorm
#' @importFrom ashr my_e2truncnorm
compute_summary_results_laplace = function(x,s,w,a){
  wpost <- wpost_laplace(x,s,w,a)
  l=lambda(x,s,a)
  posterior_mean = wpost* (l * my_etruncnorm(0,Inf,x-s^2*a,s) + (1-l)*my_etruncnorm(-Inf,0,x+s^2*a,s))
  posterior_mean2 = wpost* (l * my_e2truncnorm(0,Inf,x-s^2*a,s) + (1-l)*my_e2truncnorm(-Inf,0,x+s^2*a,s))
  return(data.frame(posterior_mean=posterior_mean,posterior_mean2=posterior_mean2))
}

#
#  Calculate the posterior weight for non-zero effect
#
#' @importFrom stats dnorm
wpost_laplace <- function(x, s, w, a)
{
  if(w==0){return(rep(0,length(x)))}
  if(w==1){return(rep(1,length(x)))}
  lg = logg_laplace(x,s,a)
  lf = dnorm(x,0,s,log=TRUE)
  return(w/(w+(1-w)*exp(lf-lg)))
}


# computes the lambda function equation (2.7) from Kan Xu's thesis, which is posterior probability of being negative
# given a non-zero effect
lambda = function(x,s,a){
  lm1 = -a*x  + pnorm(x/s - s*a,log.p = TRUE)
  lm2 =  a*x +  pnorm(x/s + s*a,log.p = TRUE,lower.tail = FALSE)
  m = pmax(lm1,lm2)
  lm1 = lm1-m
  lm2 = lm2-m
  return(exp(lm1)/(exp(lm1)+exp(lm2)))
}
