#' @title Compute mle of w and a for data (x,s) under point normal prior
#' @description Paragraph-length description goes here.
#' @details Does simple mle (currently using optim)
#' @param x observations
#' @param s standard deviations
#' @param startpar initialization
#' @param control list of parameters to be passed to optim
#' @importFrom stats optim
mle_normal <- function(x, s, startpar = NULL, control=NULL) {

  # get some reasonable limits on sigma (standard deviation of normal, or 1/sqrt(a))
  sigmamin = min(s)/10
  if (all(x^2 <= s^2)) {
    sigmamax = 8 * sigmamin
  } else {
    sigmamax = 2 * sqrt(max(x^2 - s^2))
  }

  lo  <-  c(0,1/sigmamax^2)
  hi  <-  c(1,1/sigmamin^2)
  if(is.null(startpar)){
    startpar  <- c(0.5,2/(sigmamax^2+sigmamin^2))
  }

  uu <- optim(startpar, function(par,x,s){-loglik_normal(x,s,par[1],par[2])}, method="L-BFGS-B",
                lower = lo, upper = hi, x = x, s = s, control=control)
  uu_par <- uu$par

  return(list(pi0=1-uu_par[1], a=uu_par[2], val = uu$value))
}



#do optimization of parameters on log scale
#
#' @importFrom stats optim
mle_normal_logscale <- function(x, s,startpar=NULL,control=NULL) {

  maxvar = max(x^2 - s^2) #get upper bound on variance estimate

  if(maxvar<0){ #deal with case where everything is smaller than expected (null)
    return(list(pi0=1, a=1)) # note that a is irrelevant if pi0=1
  }

  # set default starting point. This point is chosen based on the model
  # where there is a single non-null value, based on the intuition that
  # this is the case that is "hardest" to get right
  if(is.null(startpar)){
    startpar  <- c(log(1/length(x)),-log(maxvar))
  }

  uu <- optim(startpar,
              function(par,x,s){-loglik_normal(x,s,exp(par[1])/(1+exp(par[1])),
                                                exp(par[2]))},
              method="L-BFGS-B",x = x, s = s, control=control)
  uu_par <- uu$par

  return(list(pi0=1/(1+exp(uu_par[1])), a=exp(uu_par[2]), val= uu$value))
}


#do optimization of parameters on log scale
#
#' @importFrom stats optim
mle_normal_logscale_grad <- function(x, s,startpar = NULL, control=NULL) {

  maxvar = max(x^2 - s^2) #get upper bound on variance estimate

  if(maxvar<0){ #deal with case where everything is smaller than expected (null)
    return(list(pi0=1, a=1)) # note that a is irrelevant if pi0=1
  }


  # set default starting point. This point is chosen based on the model
  # where there is a single non-null value, based on the intuition that
  # this is the case that is "hardest" to get right
  if(is.null(startpar)){
    startpar  <- c(log(1/length(x)),-log(maxvar))
  }

  uu <- try(optim(startpar,
              function(par,x,s){-loglik_normal(x,s,1- (1/(1+exp(par[1]))),
                                               exp(par[2]))},
              gr= function(par,x,s){grad_negloglik_logscale_normal(x,s,1- (1/(1+exp(par[1]))),
                                                            exp(par[2]))},
              method="L-BFGS-B",x = x, s = s, control=control), silent=TRUE)

  #if optimization fails, try again with some limits; this should not really
  #happen but in preliminary testing sometimes we say optim complain of infinite
  # values, possibly because of extreme values of the parameters?
  if(class(uu)=="try-error"){
    n = length(x)
    minvar = (min(s)/10)^2

    lo  <-  c(log(1/n),-log(maxvar))
    hi  <-  c(log(n),-log(minvar))
    uu <- try(optim(startpar,
                    function(par,x,s){-loglik_normal(x,s,1- (1/(1+exp(par[1]))),
                                                     exp(par[2]))},
                    gr= function(par,x,s){grad_negloglik_logscale_normal(x,s,1- (1/(1+exp(par[1]))),
                                                                         exp(par[2]))},
                    method="L-BFGS-B",lower=lo, upper=hi, x = x, s = s, control=control))
  }
  if(class(uu)=="try-error"){
    saveRDS(list(startpar=startpar,x=x,s=s,control=control),"temp_debug.RDS")
    stop("optim failed to converge; debug information saved to temp_debug.RDS")
  }

  uu_par <- uu$par

  return(list(pi0=1/(1+exp(uu_par[1])), a=exp(uu_par[2]), val = uu$value))
}



