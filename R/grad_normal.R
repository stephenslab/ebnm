#
# Computes gradient of negloglik with respect to mu and t:=log(a)
#
#' @importFrom stats dnorm
grad_negloglik_logscale_normal = function(x, s, mu, t) {
  grad_mu = -sum((x - mu) / (exp(-t) + s^2))
  grad_t = -0.5 * exp(-t) * sum((1 / (exp(-t) + s^2)) * (1 - (((x - mu)^2) / (exp(-t) + s^2))))
  grad = c(grad_mu, grad_t)
  return(grad)
}
