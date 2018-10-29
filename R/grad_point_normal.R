#
# Computes gradient with respect to logit(w) and log(a) and mu
#
grad_negloglik_logscale_point_normal = function(x, s, w, a, mu) {
  grad = grad_negloglik_point_normal(x, s, w, a, mu)
  grad[1] = grad[1] * (w * (1 - w))
  grad[2] = grad[2] * a
  grad[3] = grad[3]
  return(grad)
}

#' @importFrom stats dnorm
#'
grad_negloglik_point_normal = function(x, s, w, a, mu) {
  s[s==0] = .Machine$double.eps # avoid numeric problems when s = 0
  l = vloglik_point_normal(x, s, w, a, mu)
  lf = dnorm(x, mu, s, log = TRUE)
  lg = dnorm(x, mu, sqrt(s^2 + 1/a), log = TRUE)

  grad_w = sum(exp(lf - l) - exp(lg - l))

  g_over_l = exp(lg - l)
  grad_a = -w * sum(g_over_l * grad_lg_point_normal(x, s, a, mu))
  
  f_over_l = exp(lf - l)
  grad_mu = -sum((w * g_over_l * (x - mu) / (s^2 + (1 / a))) + ((1 - w) * f_over_l * (x - mu) / (s^2)))

  return(c(grad_w, grad_a, grad_mu))
}

grad_lg_point_normal = function(x, s, a, mu) {
  vinv = (1 / (s^2 + 1/a))
  return(0.5 * (1 / a^2) * (vinv - (x - mu)^2 * vinv^2))
}
