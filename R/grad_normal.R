#
# Computes gradient of negloglik with respect to mu and log(a)
#
grad_negloglik_logscale_normal = function(x, s, mu, a) {
  vinv = (1 / (s^2 + 1/a))
  
  grad_mu = -sum((x - mu) * vinv)
  grad_a = 0.5 * (1/a) * sum(((x - mu) * vinv)^2 - vinv)
  grad = c(grad_mu, grad_a)
  return(grad)
}
