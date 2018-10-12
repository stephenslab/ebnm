#
# Computes gradient with respect to logit(w) and log(a)
#
grad_negloglik_logscale_point_normal = function(x, s, w, a) {
  grad = grad_negloglik_point_normal(x, s, w, a)
  grad[1] = grad[1] * (w * (1 - w))
  grad[2] = grad[2] * a
  return(grad)
}

#' @importFrom stats dnorm
#'
grad_negloglik_point_normal = function(x, s, w, a) {
  s[s==0] = .Machine$double.eps # avoid numeric problems when s = 0
  l = vloglik_point_normal(x, s, w, a)
  lf = dnorm(x, 0, s, log = TRUE)
  lg = dnorm(x, 0, sqrt(s^2 + 1/a), log = TRUE)
  
  grad_w = sum(exp(lf - l) - exp(lg - l))
  
  g_over_l = exp(lg - l)
  grad_a = -w * sum(g_over_l * grad_lg_point_normal(x, s, a))
  
  return(c(grad_w, grad_a))
}

grad_lg_point_normal = function(x, s, a) {
  vinv = (1 / (s^2 + 1/a))
  return(0.5 * (1 / a^2) * (vinv - x^2 * vinv^2))
}
