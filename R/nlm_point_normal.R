nlm_pn_fn <- function(par, fix_pi0, fix_a, fix_mu, alpha, beta, mu,
                      n0, n1, sum1, n2, x, s2, z, sum_z) {
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
    sum_z <- sum(z)
  }

  logist.alpha  <- 1 / (1 + exp(-alpha)) # scalar
  logist.nalpha <- 1 / (1 + exp(alpha))
  logist.beta   <- 1 / (1 + s2 * exp(-beta)) # scalar or vector
  logist.nbeta  <- 1 / (1 + exp(beta) / s2)

  y <- 0.5 * (z * logist.beta + log(logist.nbeta)) # vector
  logist.y  <- 1 / (1 + exp(alpha - y)) # vector
  logist.ny <- 1 / (1 + exp(y - alpha))

  # Negative log likelihood.
  C <- pmax(y, alpha)
  nllik <- -n0 * log(logist.alpha) - (n1 + n2) * (log(logist.nalpha))
  nllik <- nllik + 0.5 * (n1 * beta + sum1 * exp(-beta) + sum_z)
  nllik <- nllik - sum(log(exp(y - C) + exp(alpha - C)) + C)

  # Gradient.
  grad <- numeric(length(par))
  i <- 1
  if (!fix_pi0) {
    grad[i] <- -n0 * logist.nalpha + (n1 + n2) * logist.alpha - sum(logist.ny)
    i <- i + 1
  }
  if (!fix_a) {
    dy.dbeta <- 0.5 * (z * logist.beta * logist.nbeta - logist.beta)
    grad[i] <- 0.5 * (n1 - sum1 * exp(-beta)) - sum(logist.y * dy.dbeta)
    i <- i + 1
  }
  if (!fix_mu) {
    dy.dmu <- (mu - x) * logist.beta / s2
    grad[i] <- sum((mu - x) / s2) - sum(logist.y * dy.dmu)
  }

  # Hessian.
  hess <- matrix(nrow = length(par), ncol = length(par))
  i <- 1
  if (!fix_pi0) {
    tmp <- (n0 + n1 + n2) * logist.alpha * logist.nalpha
    hess[i, i] <- tmp - sum(logist.y * logist.ny)
    j <- i + 1
    if (!fix_a) {
      hess[i, j] <- hess[j, i] <- sum(logist.y * logist.ny * dy.dbeta)
      j <- j + 1
    }
    if (!fix_mu) {
      hess[i, j] <- hess[j, i] <- sum(logist.y * logist.ny * dy.dmu)
    }
    i <- i + 1
  }
  if (!fix_a) {
    d2y.dbeta2 <- 0.5 * ((z * (logist.nbeta - logist.beta) - 1)
                         * logist.beta * logist.nbeta)
    tmp <- 0.5 * sum1 * exp(-beta) - sum(logist.y * logist.ny * dy.dbeta^2)
    hess[i, i] <- tmp - sum(logist.y * d2y.dbeta2)
    j <- i + 1
    if (!fix_mu) {
      d2y.dbetadmu <- (mu - x) * logist.beta * logist.nbeta / s2
      tmp <- -sum(logist.y * logist.ny * dy.dbeta * dy.dmu)
      hess[i, j] <- hess[j, i] <- tmp - sum(logist.y * d2y.dbetadmu)
    }
    i <- i + 1
  }
  if (!fix_mu) {
    tmp <- sum(1 / s2) - sum(logist.y * logist.ny * dy.dmu^2)
    hess[i, i] <- tmp - sum(logist.y * logist.beta / s2)
  }

  attr(nllik, "gradient") <- grad
  attr(nllik, "hessian") <- hess
  return(nllik)
}

pn_llik_from_optval <- function(optval, n1, n2, s2) {
  if (length(s2) == 1) {
    sum.log.s2 <- n2 * log(s2)
  } else {
    sum.log.s2 <- sum(log(s2))
  }
  return(-optval - 0.5 * ((n1 + n2) * log(2 * pi) + sum.log.s2))
}
