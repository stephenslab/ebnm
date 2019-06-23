# Computes MLE for g under point-normal prior.
#
#' @importFrom stats optim
#'
mle_point_laplace <- function(x, s, startpar = NULL, control = NULL) {
  startpar <- pl_startpar(x, s)

  optres <- try(optim(startpar, pl_fn, pl_gr,
                      x = x, s = s,
                      method = "L-BFGS-B", control = control),
                silent = TRUE)

  if (inherits(optres, "try-error") || optres$convergence != 0) {
    warning("First optimization attempt failed. Retrying with bounds.")
    hilo <- pl_hilo(x, s)
    optres <- optim(startpar, pl_fn, pl_gr,
                    x = x, s = s,
                    method = "L-BFGS-B", control = control,
                    lower = hilo$lo, upper = hilo$hi)
  }

  retlist <- pl_g_from_optpar(optres$par)
  retlist$val <- optres$val

  return(retlist)
}

# Initial values.
pl_startpar <- function(x, s) {
  # Defaults for logit(pi0) and log(a).
  return(c(0, -log(mean(x^2))))
}

# Negative log likelihood.
pl_fn <- function(par, x, s) {
  return(-loglik_point_laplace(x, s, exp(par[1]) / (1 + exp(par[1])), exp(par[2])))
}

# Gradient of the negative log likelihood.
pl_gr <- function(par, x, s) {
  w <- exp(par[1]) / (1 + exp(par[1]))
  a <- exp(par[2])

  l = vloglik_point_laplace(x, s, w, a)
  lf = dnorm(x/s, log = TRUE)
  lg = logg_laplace(x, s, a)
  grad_w = sum(exp(lf - l) - exp(lg - l))

  f_over_g = exp(lf - lg)

  lg1 <- -a * x + pnorm((x - s^2 * a) / s, log.p = TRUE)
  lg2 <-  a * x + pnorm((x + s^2 * a) / s, log.p = TRUE, lower.tail = FALSE)
  grad_lg1 <- -x - s * exp(dnorm(x/s - s * a, log = TRUE)
                           - pnorm(x/s - s * a, log.p = TRUE))
  grad_lg2 <-  x - s * exp(dnorm(x/s + s * a, log = TRUE)
                           - pnorm(x/s + s * a, log.p = TRUE, lower.tail = FALSE))
  weight <- 1 / (1 + exp(lg2 - lg1))
  grad_lg <- 1 / a + a * s^2 + weight * grad_lg1 + (1 - weight) * grad_lg2
  grad_a <- -w * sum(grad_lg / ((1 - w) * f_over_g + w))

  grad <- c(grad_w, grad_a)
  grad[1] = grad[1] * (w * (1 - w))
  grad[2] = grad[2] * a

  return(grad)
}

pl_g_from_optpar <- function(par) {
  return(list(pi0 = 1 / (1 + exp(par[1])), a = exp(par[2])))
}

# Upper and lower bounds for optim in case the first attempt fails.
pl_hilo <- function(x, s) {
  n <- length(x)
  maxvar <- max(x^2 - s^2)
  minvar <- (min(s) / 10)^2

  lo <- c(log(1/n), -log(maxvar))
  hi <- c(log(n), -log(minvar))
  return(list(lo = lo, hi = hi))
}
