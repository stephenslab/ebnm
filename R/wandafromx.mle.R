#' @title Compute mle of w and a for data (x,s)
#' @details Does simple mle (currently using optim)
#' Compared with the original wandafromx this function dispenses with the
#' complication of the thresholds
wandafromx.mle <- function(x, s) {

  # get some reasonable limits on sigma (standard deviation of laplace, or 1/sqrt(a))
  sigmamin = min(s)/10
  if (all(x^2 <= s^2)) {
    sigmamax = 8 * sigmamin
  } else {
    sigmamax = 2 * sqrt(max(x^2 - s^2))
  }

  lo  <-  c(0,1/sigmamax^2)
  hi  <-  c(1,1/sigmamin^2)
  startpar  <- c(0.5,2/(sigmamax^2+sigmamin^2))

  uu <- optim(startpar, function(par,x,s){-loglik.laplace(x,s,par[1],par[2])}, method="L-BFGS-B",
                lower = lo, upper = hi, x = x, s = s)
  uu <- uu$par

  return(list(w=uu[1], a=uu[2]))
}

wandafromx.mle.grad <- function(x, s) {

  # get some reasonable limits on sigma (standard deviation of laplace, or 1/sqrt(a))
  sigmamin = min(s)/10
  if (all(x^2 <= s^2)) {
    sigmamax = 8 * sigmamin
  } else {
    sigmamax = 2 * sqrt(max(x^2 - s^2))
  }

  lo  <-  c(0,1/sigmamax^2)
  hi  <-  c(1,1/sigmamin^2)
  startpar  <- c(0.5,2/(sigmamax^2+sigmamin^2))

  uu <- optim(startpar, function(par,x,s){-loglik.laplace(x,s,par[1],par[2])},
              gr= function(par,x,s){grad_negloglik.laplace(x,s,par[1],par[2])}, method="L-BFGS-B",
              lower = lo, upper = hi, x = x, s = s)
  uu <- uu$par

  return(list(w=uu[1], a=uu[2]))
}





