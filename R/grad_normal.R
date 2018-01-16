# gradient of negative likelihood with respect to w=(1-pi0)
#
#' @importFrom stats dnorm
grad_negloglik_w_normal = function(x,s,w,a){
  l = vloglik_normal(x,s,w,a)
  lf = dnorm(x,0,s,log=TRUE)
  lg = dnorm(x,0,sqrt(s^2+1/a),log=TRUE)
  sum(exp(lf-l)-exp(lg-l))
}

#' @importFrom stats dnorm
grad_negloglik_a_normal = function(x,s,w,a){
  lg = dnorm(x,0,sqrt(s^2 + 1/a),log=TRUE)
  l = vloglik_normal(x,s,w,a)
  g_over_l = exp(lg-l)
  -w*sum(g_over_l*grad_lg_normal(x,s,a))
}

#check with
#numDeriv::grad(function(a){-loglik_normal(x,s,w,a)},a)
#numDeriv::grad(function(w){-loglik_normal(x,s,w,a)},w)

# Combines the two above.
#
#' @importFrom stats dnorm
grad_negloglik_normal  = function(x,s,w,a){
  l = vloglik_normal(x,s,w,a)
  lf = dnorm(x,0,s,log=TRUE)
  lg = dnorm(x,0,sqrt(s^2+1/a),log=TRUE)
  grad_w = sum(exp(lf-l)-exp(lg-l))

  g_over_l = exp(lg-l)
  grad_a = -w*sum(g_over_l*grad_lg_normal(x,s,a))

  c(grad_w,grad_a)
}



# computes gradient with respect to logit(w) and log(a)
grad_negloglik_logscale_normal  = function(x,s,w,a){
  grad = grad_negloglik_normal(x,s,w,a)
  grad[1] = grad[1] * (w*(1-w))
  grad[2] = grad[2] * a
  return(grad)
}

# set.seed(1)
# x = rnorm(100)
# s = rgamma(100,1,1)
# w=0.7
# a=4.2
# grad_negloglik_normal(x,s,w,a)
# #[1] -35.20593
# numDeriv::grad(function(w){-loglik_normal(x,s,w,a)},w)
# #[1] -35.20593
#
#
# grad_negloglik(x,s,w,a)
# numDeriv::grad(function(w){-loglik_normal(x,s,w,a)},w)
# numDeriv::grad(function(a){-loglik_normal(x,s,w,a)},a)
#

#' @importFrom stats dnorm
lg_normal = function(x,s,a) {
  dnorm(x,0,sqrt(s^2+1/a),log=TRUE)
}

grad_lg_normal = function(x,s,a){
  vinv = (1/(s^2 + 1/a))
  0.5 * (1/a^2) * (vinv - x^2*vinv^2)
}

# could be useful...
#To save myself some work, I am calling the methods through the package optimrx
#which allows multiple methods to be applied via the single function opm()
