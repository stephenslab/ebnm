#' @title Compute mle of w and a for data (x,s) under point normal prior
#' @details Does simple mle (currently using optim)
mle_normal <- function(x, s,control=NULL) {

  # get some reasonable limits on sigma (standard deviation of normal, or 1/sqrt(a))
  sigmamin = min(s)/10
  if (all(x^2 <= s^2)) {
    sigmamax = 8 * sigmamin
  } else {
    sigmamax = 2 * sqrt(max(x^2 - s^2))
  }

  lo  <-  c(0,1/sigmamax^2)
  hi  <-  c(1,1/sigmamin^2)
  startpar  <- c(0.5,2/(sigmamax^2+sigmamin^2))

  uu <- optim(startpar, function(par,x,s){-loglik_normal(x,s,par[1],par[2])}, method="L-BFGS-B",
                lower = lo, upper = hi, x = x, s = s, control=control)
  uu <- uu$par

  return(list(pi0=1-uu[1], a=uu[2]))
}



#do optimization of parameters on log scale
mle_normal_logscale <- function(x, s,control=NULL) {

  # get some reasonable limits on sigma (standard deviation of laplace, or 1/sqrt(a))
  sigmamin = min(s)/10
  if (all(x^2 <= s^2)) {
    sigmamax = 8 * sigmamin
  } else {
    sigmamax = 2 * sqrt(max(x^2 - s^2))
  }

  startpar  <- c(0,log(2/(sigmamax^2+sigmamin^2)))

  uu <- optim(startpar,
              function(par,x,s){-loglik_normal(x,s,exp(par[1])/(1+exp(par[1])),
                                                exp(par[2]))},
              method="L-BFGS-B",x = x, s = s, control=control)
  uu <- uu$par

  return(list(pi0=1/(1+exp(uu[1])), a=exp(uu[2])))
}


#do optimization of parameters on log scale
mle_normal_logscale_grad <- function(x, s,control=NULL) {

  # get some reasonable limits on sigma (standard deviation of laplace, or 1/sqrt(a))
  sigmamin = min(s)/10
  if (all(x^2 <= s^2)) {
    sigmamax = 8 * sigmamin
  } else {
    sigmamax = 2 * sqrt(max(x^2 - s^2))
  }

  startpar  <- c(0,log(2/(sigmamax^2+sigmamin^2)))

  uu <- optim(startpar,
              function(par,x,s){-loglik_normal(x,s,exp(par[1])/(1+exp(par[1])),
                                               exp(par[2]))},
              gr= function(par,x,s){grad_negloglik_logscale_normal(x,s,exp(par[1])/(1+exp(par[1])),
                                                            exp(par[2]))},
              method="L-BFGS-B",x = x, s = s, control=control)
  uu <- uu$par

  return(list(pi0=1/(1+exp(uu[1])), a=exp(uu[2])))
}
