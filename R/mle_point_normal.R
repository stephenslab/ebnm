# Computes MLE for g under point-normal prior.
mle_point_normal <- function(x, s, g, control, fix_pi0, fix_a, fix_mu) {
  startpar <- pn_startpar(x, s, g, fix_pi0, fix_a, fix_mu)

  if (fix_pi0) {
    alpha <- -log(1 / g$pi0 - 1)
  } else {
    alpha <- NULL
  }
  if (fix_a) {
    beta <- -log(g$a)
  } else {
    beta <- NULL
  }
  if (fix_mu) {
    mu <- g$mu
  } else {
    mu <- NULL
  }

  if (any(s == 0)) {
    which.s0 <- which(s == 0)
    which.x.nz <- which(x[which.s0] != mu)
    n0 <- length(which.s0) - length(which.x.nz)
    n1 <- length(which.x.nz)
    sum1 <- sum((x[which.s0[which.x.nz]] - mu)^2)
    x <- x[-which.s0]
    s <- s[-which.s0]
  } else {
    n0 <- 0
    n1 <- 0
    sum1 <- 0
  }
  n2 <- length(x)

  s2 <- s^2

  if (fix_mu) {
    z <- (x - mu)^2 / s2
    sum.z <- sum(z)
  } else {
    z <- NULL
    sum.z <- NULL
  }

  optres <- try(optim(startpar, pn_fn, pn_gr,
                      fix_pi0 = fix_pi0, fix_a = fix_a, fix_mu = fix_mu,
                      alpha = alpha, beta = beta, mu = mu,
                      n0 = n0, n1 = n1, sum1 = sum1, n2 = n2,
                      x = x, s2 = s2, z = z, sum.z = sum.z,
                      method = "L-BFGS-B", control = control),
                silent = TRUE)

  # TODO: is this necessary?
  if (inherits(optres, "try-error") || optres$convergence != 0) {
    warning("First optimization attempt failed. Retrying with bounds.")
    hilo <- pn_hilo(x, s, fix_pi0, fix_a, fix_mu)
    optres <- optim(startpar, pn_fn, pn_gr,
                    fix_pi0 = fix_pi0, fix_a = fix_a, fix_mu = fix_mu,
                    alpha = alpha, beta = beta, mu = mu,
                    n0 = n0, n1 = n1, sum1 = sum1, n2 = n2,
                    x = x, s2 = s2, z = z, sum.z = sum.z,
                    method = "L-BFGS-B", control = control,
                    lower = hilo$lo, upper = hilo$hi)
  }

  retlist <- pn_g_from_optpar(optres$par, g, fix_pi0, fix_a, fix_mu)
  retlist$val <- pn_llik_from_optval(optres$value, n1, n2, s2)

  return(retlist)
}

# Initial values.
pn_startpar <- function(x, s, g, fix_pi0, fix_a, fix_mu) {
  startpar <- numeric(0)

  if (!fix_pi0) {
    if (!is.null(g$pi0) && g$pi0 > 0 && g$pi0 < 1) {
      startpar <- c(startpar, -log(1 / g$pi0 - 1))
    } else {
      startpar <- c(startpar, 0) # default for logit(pi0)
    }
  }

  if (!fix_a) {
    if (!is.null(g$a)) {
      startpar <- c(startpar, -log(g$a))
    } else {
      startpar <- c(startpar, log(mean(x^2))) # default for -log(a)
    }
  }

  if (!fix_mu) {
    if (!is.null(g$mu)) {
      startpar <- c(startpar, g$mu)
    } else {
      startpar <- c(startpar, mean(x)) # default for mu
    }
  }

  return(startpar)
}

# Negative log likelihood.
pn_fn <- function(par, fix_pi0, fix_a, fix_mu, alpha, beta, mu,
                  n0, n1, sum1, n2, x, s2, z, sum.z) {
  i <- 1
  if (!fix_pi0) {
    alpha <- par[i]
    i <- i + 1
  }
  if (!fix_a) {
    beta <- par[i]
    i <- i + 1
  }
  if (!fix_mu) {
    mu <- par[i]
    z <- (x - mu)^2 / s2
    sum.z <- sum(z)
  }

  y <- (z / (1 + s2 * exp(-beta)) - log(1 + exp(beta) / s2)) / 2
  C <- pmax(y, alpha)

  nllik <- n0 * log(1 + exp(-alpha)) + (n1 + n2) * (log(1 + exp(alpha)))
  nllik <- nllik + n1 * beta / 2 + sum1 * exp(-beta) / 2 + sum.z / 2
  nllik <- nllik - sum(log(exp(y - C) + exp(alpha - C)) + C)

  return(nllik)
}

# Gradient of the negative log likelihood.
pn_gr <- function(par, fix_pi0, fix_a, fix_mu, alpha, beta, mu,
                  n0, n1, sum1, n2, x, s2, z, sum.z) {
  i <- 1
  if (!fix_pi0) {
    alpha <- par[i]
    i <- i + 1
  }
  if (!fix_a) {
    beta <- par[i]
    i <- i + 1
  }
  if (!fix_mu) {
    mu <- par[i]
    z <- (x - mu)^2 / s2
  }

  tmp1 <- 1 / (1 + s2 * exp(-beta))
  tmp2 <- 1 / (1 + exp(beta) / s2)
  y <- (z * tmp1 + log(tmp2)) / 2

  grad <- numeric(0)
  if (!fix_pi0) {
    tmp <- -n0 / (1 + exp(alpha)) + (n1 + n2) / (1 + exp(-alpha))
    grad <- c(grad, tmp - sum(1 / (1 + exp(y - alpha))))
  }

  if (!fix_a) {
    tmp <- n1 / 2 - exp(-beta) * sum1 / 2
    grad <- c(grad, tmp - sum((z * tmp1 * tmp2 - tmp1) / (1 + exp(alpha - y))) / 2)
  }

  if (!fix_mu) {
    tmp <- sum((mu - x) / s2)
    grad <- c(grad, tmp + sum(tmp1 * (x - mu) / (1 + exp(alpha - y)) / s2))
  }

  return(grad)
}

# Pull pi0, a, and mu out of the optimization results.
pn_g_from_optpar <- function(optpar, g, fix_pi0, fix_a, fix_mu) {
  opt_g <- list()

  i <- 1
  if (fix_pi0) {
    opt_g$pi0 <- g$pi0
  } else {
    opt_g$pi0 <- 1 / (exp(-optpar[i]) + 1)
    i <- i + 1
  }

  if (fix_a) {
    opt_g$a <- g$a
  } else {
    opt_g$a <- exp(-optpar[i])
    i <- i + 1
  }

  if (fix_mu) {
    opt_g$mu <- g$mu
  } else {
    opt_g$mu <- optpar[i]
  }

  return(opt_g)
}

pn_llik_from_optval <- function(optval, n1, n2, s2) {
  return(-optval - 0.5 * ((n1 + n2) * log(2 * pi) + sum(log(s2))))
}

# Upper and lower bounds for optim in case the first attempt fails.
pn_hilo <- function(x, s, fix_pi0, fix_a, fix_mu) {
  lo <- numeric(0)
  hi <- numeric(0)

  if (!fix_pi0) {
    n <- length(x)
    lo <- c(lo, log(1 / n))
    hi <- c(hi, log(n))
  }

  if (!fix_a) {
    maxvar <- (max(x) - min(x))^2
    minvar <- (min(s) / 10)^2
    if (minvar < 1e-8) {
      minvar <- 1e-8
    }
    lo <- c(lo, log(minvar))
    hi <- c(hi, log(maxvar))
  }

  if (!fix_mu) {
    lo <- c(lo, min(x) - 3 * max(s) - 3 * max(abs(x)))
    hi <- c(hi, max(x) + 3 * max(s) + 3 * max(abs(x)))
  }

  return(list(lo = lo, hi = hi))
}
