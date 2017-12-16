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



# gradient of negative likelihood with respect to w=(1-pi0)
grad_w = function(x,s,w,a){
  l = vloglik.laplace(x,s,w,a)
  lf = dnorm(x,0,s,log=TRUE)
  lg = logg.laplace(x,s,a)
  lfac = pmax(lf,lg)
  lnum = lfac + log(exp(lf-lfac)-exp(lg-lfac)) #numerator
  ldenom = l
  sum(exp(lnum-ldenom))
}

# Testing:
# negloglik.w = function(w){
#   -loglik.laplace(x,s,w,a)
# }
# x=1
# s=1
# a=0.5
# > numDeriv::grad(negloglik.w,0.5)
# [1] 0.4691994
# > grad_w(x,s,0.5,a)
# [1] 0.4691994


# gradient of g function
grad_g = function(x,s,w,a){
  x = abs(x) # allows us to do all calculations as if x is positive
  xma = (x - s^2*a)/s
  xpa = (x + s^2*a)/s

  lterm1 = log(s) - a*x + dnorm(xma,log=TRUE)
  lterm2 = log(x) - a*x + pnorm(xma,log=TRUE)
  lfac = pmax(lterm1,lterm2)
  logA = lfac + log(exp(lterm1-lfac)+exp(lterm2-lfac))

  lterm1 = log(s) + a*x + dnorm(xpa,log=TRUE)
  lterm2 = log(x) + a*x + pnorm(xpa,log=TRUE,lower.tail=FALSE)
  lfac = pmax(lterm1,lterm2)
  B = exp(lfac) * (exp(lterm2-lfac)- exp(lterm1-lfac))

  C = B - exp(logA)

  return((1/a + a*s^2)*exp(logg.laplace(x,s,a)) + 0.5*a*exp(0.5*a^2 * s^2)*C)
}

g.laplace = function(a){
   logg.laplace(x,s,a)
}
x=1
s=1
w=0.5
numDeriv::grad(g.laplace,0.5)
# [1] 0.4691994
# > grad_w(x,s,0.5,a)
# [1] 0.4691994

lg1 = function(x,s,a){ -a*x + pnorm((x-s^2*a)/s,log=TRUE)}
grad_lg1 = function(x,s,a){
  -x -s*dnorm(x/s-s*a)/pnorm(x/s-s*a)
}

lg2 = function(x,s,a){a*x + pnorm((x+s^2*a)/s,lower.tail = FALSE,log=TRUE)}
grad_lg2= function(x,s,a){
  x - s*dnorm(x/s+s*a)/pnorm(x/s+s*a,lower.tail=FALSE)
}

grad_lg = function(x,s,a){
  1/a + a*s^2 +
    (grad_lg1(x,s,a)*exp(lg1(x,s,a)) + grad_lg2(x,s,a)*exp(lg2(x,s,a))) /
    (exp(lg1(x,s,a))+exp(lg2(x,s,a)))
}

grad_lg(x,s,0.5)
numDeriv::grad(function(a){logg.laplace(x,s,a)},0.5)

#grad_lg2(x,s,0.5)
#numDeriv::grad(lg2,0.5)


# test
#> x
# [1] -1
# > a
# [1] 0.5
# > w
# [1] 0.5
#> grad_lg1(0.5)
#[1] -0.9386772
#> numDeriv::grad(lg1,0.5)
#[1] -0.9386772


# could be useful...
#To save myself some work, I am calling the methods through the package optimrx
#which allows multiple methods to be applied via the single function opm()
