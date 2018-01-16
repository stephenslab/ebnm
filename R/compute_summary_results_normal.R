compute_summary_results_normal = function(x,s,w,a){
  wpost <- wpost_normal(x,s,w,a)
  pmean_cond = x/(1+s^2*a) # posterior mean from standard normal analysis
  pvar_cond = s^2/(1+s^2*a) # posterior variance from standard normal analysis
  PosteriorMean = wpost* pmean_cond
  PosteriorMean2 = wpost* (pmean_cond^2 + pvar_cond)
  return(data.frame(PosteriorMean=PosteriorMean,PosteriorMean2=PosteriorMean2))
}

#
#  Calculate the posterior weight for non-zero effect
#
#' @importFrom stats dnorm
wpost_normal <- function(x, s, w, a)
{
  if(w==0){return(rep(0,length(x)))}
  if(w==1){return(rep(1,length(x)))}
  lg = dnorm(x,0,s+sqrt(1/a),log=TRUE)
  lf = dnorm(x,0,s,log=TRUE)
  return(w/(w+(1-w)*exp(lf-lg)))
}

