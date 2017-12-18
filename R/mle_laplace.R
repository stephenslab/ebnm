#' @title Compute mle of w and a for data (x,s) under point laplace prior
#' @details Does simple mle (currently using optim)
#' @param x observations
#' @param s standard deviations
#' @param startpar initialization
#' @param control list of parameters to be passed to optim
#' Compared with the original wandafromx this function dispenses with the
#' complication of the thresholds
mle_laplace <- function(x, s,startpar=NULL,control=NULL) {

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

  uu <- optim(startpar, function(par,x,s){-loglik_laplace(x,s,par[1],par[2])}, method="L-BFGS-B",
                lower = lo, upper = hi, x = x, s = s, control=control)
  uu <- uu$par

  return(list(pi0=1-uu[1], a=uu[2]))
}

mle_laplace_grad <- function(x, s, startpar=NULL,control=NULL) {

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
  uu <- optim(startpar, function(par,x,s){-loglik_laplace(x,s,par[1],par[2])},
              gr= function(par,x,s){grad_negloglik(x,s,par[1],par[2])}, method="L-BFGS-B",
              lower = lo, upper = hi, x = x, s = s, control=control)
  uu <- uu$par

  return(list(pi0=1-uu[1], a=uu[2]))
}

#do optimization of parameters on log scale
mle_laplace_logscale <- function(x, s,startpar=NULL,control=NULL) {

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
  n = length(x)
  minvar = (min(s)/10)^2
  lo  <-  c(log(1/n),-log(maxvar))
  hi  <-  c(log(n),-log(minvar))

  uu <- optim(startpar,
              function(par,x,s){-loglik_laplace(x,s,exp(par[1])/(1+exp(par[1])),
                                                exp(par[2]))},
              method="L-BFGS-B",
              lower = lo, upper = hi, x = x, s = s, control=control)
  uu <- uu$par

  return(list(pi0=1/(1+exp(uu[1])), a=exp(uu[2])))
}


#do optimization of parameters on log scale
mle_laplace_logscale_grad <- function(x, s,startpar=NULL,control=NULL) {

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
  n = length(x)
  minvar = (min(s)/10)^2
  lo  <-  c(log(1/n),-log(maxvar))
  hi  <-  c(log(n),-log(minvar))


  uu <- optim(startpar,
              function(par,x,s){-loglik_laplace(x,s,exp(par[1])/(1+exp(par[1])),
                                                exp(par[2]))},
              gr= function(par,x,s){grad_negloglik_logscale_laplace(x,s,exp(par[1])/(1+exp(par[1])),
                                                            exp(par[2]))},
              method="L-BFGS-B",
              lower = lo, upper = hi, x = x, s = s, control=control)
  uu <- uu$par

  return(list(pi0=1/(1+exp(uu[1])), a=exp(uu[2])))
}


# x = rnorm(10000) + c(rep(0,5000),rexp(5000))
# s = 1
# microbenchmark::microbenchmark(wandafromx.mle_grad(x,s),times=10)
# microbenchmark::microbenchmark(wandafromx.mle(x,s),times=10)
# microbenchmark::microbenchmark(wandafromx.mle_logscale_grad(x,s),times=10)
