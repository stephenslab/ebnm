# gradient of negative likelihood with respect to w=(1-pi0)
grad_negloglik_w = function(x,s,w,a){
  l = vloglik.laplace(x,s,w,a)
  lf = dnorm(x,0,s,log=TRUE)
  lg = logg.laplace(x,s,a)
  sum(exp(lf-l)-exp(lg-l))
}

grad_negloglik_a = function(x,s,w,a){
  lf = dnorm(x,0,s,log=TRUE)
  lg = logg.laplace(x,s,a)
  f_over_g = exp(lf-lg)
  -w * sum(grad_lg(x,s,a)/((1-w)*f_over_g+w))
}

# combines the two above
grad_negloglik  = function(x,s,w,a){
  l = vloglik.laplace(x,s,w,a)
  lf = dnorm(x,0,s,log=TRUE)
  lg = logg.laplace(x,s,a)
  grad_w = sum(exp(lf-l)-exp(lg-l))

  f_over_g = exp(lf-lg)
  grad_a = -w * sum(grad_lg(x,s,a)/((1-w)*f_over_g+w))

  c(grad_w,grad_a)
}



# computes gradient with respect to logit(w) and log(a)
grad_negloglik.logscale  = function(x,s,w,a){
  grad = grad_negloglik(x,s,w,a)
  grad[1] = grad[1] * (w*(1-w))
  grad[2] = grad[2] * a
  return(grad)
}

# set.seed(1)
# x = rnorm(100)
# s = rgamma(100,1,1)
# w=0.7
# a=4.2
# grad_negloglik_w(x,s,w,a)
# #[1] -35.20593
# numDeriv::grad(function(w){-loglik.laplace(x,s,w,a)},w)
# #[1] -35.20593
#
#
# grad_negloglik(x,s,w,a)
# numDeriv::grad(function(w){-loglik.laplace(x,s,w,a)},w)
# numDeriv::grad(function(a){-loglik.laplace(x,s,w,a)},a)
#

lg1 = function(x,s,a){ -a*x + pnorm((x-s^2*a)/s,log=TRUE)}
grad_lg1 = function(x,s,a){
  -x -s*exp(dnorm(x/s-s*a,log=TRUE)-pnorm(x/s-s*a,log=TRUE))
}

lg2 = function(x,s,a){a*x + pnorm((x+s^2*a)/s,lower.tail = FALSE,log=TRUE)}
grad_lg2= function(x,s,a){
  x - s*exp(dnorm(x/s+s*a,log=TRUE)-pnorm(x/s+s*a,lower.tail=FALSE,log=TRUE))
}

grad_lg = function(x,s,a){
  weight = 1/(1+exp(lg2(x,s,a)-lg1(x,s,a)))
  1/a + a*s^2 + weight*grad_lg1(x,s,a) + (1-weight)*grad_lg2(x,s,a)
}


grad_g = function(x,s,a){
  g = exp(logg.laplace(x,s,a))
  return(g*(1/a + a*s^2) + 0.5*a*exp(0.5*a^2*s^2)*
           (grad_lg1(x,s,a)*exp(lg1(x,s,a)) + grad_lg2(x,s,a)*exp(lg2(x,s,a))))
}


# could be useful...
#To save myself some work, I am calling the methods through the package optimrx
#which allows multiple methods to be applied via the single function opm()
