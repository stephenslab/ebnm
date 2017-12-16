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
