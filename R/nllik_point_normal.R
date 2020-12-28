# Point-normal uses the ashr class normalmix.
pn_checkg <- function(g_init, fix_g, pointmass, mode, scale, call) {
  return(check_g_init(g_init, fix_g, pointmass, mode, scale, call,
                      class_name = "normalmix", scale_name = "sd"))
}

# Point-normal parameters are logit(pi0), log(s2), and mu.
pn_gtopar <- function(g_init, fix_par) {
  fix_pi0 <- fix_par[1]
  fix_s2  <- fix_par[2]
  fix_mu  <- fix_par[3]

  if (!is.null(g_init) && length(g_init$pi) == 1) {
    par <- list(pi0 = 0,
                a = 1 / g_init$sd^2,
                mu = g_init$mean)
  } else if (!is.null(g_init) && length(g_init$pi) == 2) {
    par <- list(pi0 = g_init$pi[1],
                a = 1 / g_init$sd[2]^2,
                mu = g_init$mean[1])
  } else {
    par <- list()
    if (fix_pi0) {
      par$pi0 <- 0
    }
    if (fix_s2) {
      par$a <- 1 / scale^2
    }
    if (fix_mu) {
      par$mu <- mode
    }
  }

  return(par)
}

pn_startpar <- function(x, s, g, fix_par) {
  fix_pi0 <- fix_par[1]
  fix_s2  <- fix_par[2]
  fix_mu  <- fix_par[3]

  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }

  startpar <- numeric(0)

  if (!fix_pi0) {
    if (!is.null(g$pi0) && g$pi0 > 0 && g$pi0 < 1) {
      startpar <- c(startpar, -log(1 / g$pi0 - 1))
    } else {
      startpar <- c(startpar, 0) # default for logit(pi0)
    }
  }

  if (!fix_s2) {
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

pn_precomp <- function(x, s, g, fix_par) {
  fix_pi0 <- fix_par[1]
  fix_s2  <- fix_par[2]
  fix_mu  <- fix_par[3]

  if (fix_pi0) {
    alpha <- -log(1 / g$pi0 - 1)
  } else {
    alpha <- NULL
  }
  if (fix_s2) {
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
    which_s0 <- which(s == 0)
    which_x_nz <- which(x[which_s0] != g$mu)
    n0 <- length(which_s0) - length(which_x_nz)
    n1 <- length(which_x_nz)
    sum1 <- sum((x[which_s0[which_x_nz]] - g$mu)^2)
    x <- x[-which_s0]
    s <- s[-which_s0]
  } else {
    n0 <- 0
    n1 <- 0
    sum1 <- 0
  }
  n2 <- length(x)

  s2 <- s^2

  if (fix_mu) {
    z <- (x - g$mu)^2 / s2
    sum_z <- sum(z)
  } else {
    z <- NULL
    sum_z <- NULL
  }

  return(list(alpha = alpha, beta = beta, mu = mu, n0 = n0, n1 = n1,
              sum1 = sum1, n2 = n2, s2 = s2, z = z, sum_z = sum_z))
}

pn_nllik <- function(par, x, s, g, fix_par,
                     alpha, beta, mu, n0, n1, sum1, n2, s2, z, sum_z,
                     calc_grad, calc_hess) {
  fix_pi0 <- fix_par[1]
  fix_s2  <- fix_par[2]
  fix_mu  <- fix_par[3]

  i <- 1
  if (!fix_pi0) {
    alpha <- par[i]
    i <- i + 1
  }
  if (!fix_s2) {
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

  # Negative log likelihood.
  C <- pmax(y, alpha)
  nllik <- -n0 * log(logist.alpha) - (n1 + n2) * (log(logist.nalpha))
  nllik <- nllik + 0.5 * (n1 * beta + sum1 * exp(-beta) + sum_z)
  nllik <- nllik - sum(log(exp(y - C) + exp(alpha - C)) + C)

  if (calc_grad || calc_hess) {
    dlogist.beta  <- logist.beta * logist.nbeta

    logist.y  <- 1 / (1 + exp(alpha - y)) # vector
    logist.ny <- 1 / (1 + exp(y - alpha))

    # Gradient.
    grad <- numeric(length(par))
    i <- 1
    if (!fix_pi0) {
      grad[i] <- -n0 * logist.nalpha + (n1 + n2) * logist.alpha - sum(logist.ny)
      i <- i + 1
    }
    if (!fix_s2) {
      dy.dbeta <- 0.5 * (z * dlogist.beta - logist.beta)
      grad[i] <- 0.5 * (n1 - sum1 * exp(-beta)) - sum(logist.y * dy.dbeta)
      i <- i + 1
    }
    if (!fix_mu) {
      dy.dmu <- (mu - x) * logist.beta / s2
      grad[i] <- sum((mu - x) / s2) - sum(logist.y * dy.dmu)
    }
    attr(nllik, "gradient") <- grad
  }

  if (calc_hess) {
    dlogist.alpha <- logist.alpha * logist.nalpha
    dlogist.y <- logist.y * logist.ny

    # Hessian.
    hess <- matrix(nrow = length(par), ncol = length(par))
    i <- 1
    if (!fix_pi0) {
      tmp <- (n0 + n1 + n2) * dlogist.alpha
      hess[i, i] <- tmp - sum(dlogist.y)
      j <- i + 1
      if (!fix_s2) {
        hess[i, j] <- hess[j, i] <- sum(dlogist.y * dy.dbeta)
        j <- j + 1
      }
      if (!fix_mu) {
        hess[i, j] <- hess[j, i] <- sum(dlogist.y * dy.dmu)
      }
      i <- i + 1
    }
    if (!fix_s2) {
      d2y.dbeta2 <- 0.5 * ((z * (logist.nbeta - logist.beta) - 1) * dlogist.beta)
      tmp <- 0.5 * sum1 * exp(-beta) - sum(dlogist.y * dy.dbeta^2)
      hess[i, i] <- tmp - sum(logist.y * d2y.dbeta2)
      j <- i + 1
      if (!fix_mu) {
        d2y.dbetadmu <- (mu - x) * dlogist.beta / s2
        tmp <- -sum(dlogist.y * dy.dbeta * dy.dmu)
        hess[i, j] <- hess[j, i] <- tmp - sum(logist.y * d2y.dbetadmu)
      }
      i <- i + 1
    }
    if (!fix_mu) {
      tmp <- sum(1 / s2 - dlogist.y * dy.dmu^2)
      hess[i, i] <- tmp - sum(logist.y * logist.beta / s2)
    }
    attr(nllik, "hessian") <- hess
  }

  return(nllik)
}

pn_gfromopt <- function(optpar, optval, g, fix_par,
                        alpha, beta, mu, n0, n1, sum1, n2, s2, z, sum_z) {
  fix_pi0 <- fix_par[1]
  fix_s2   <- fix_par[2]
  fix_mu  <- fix_par[3]

  opt_g <- list()

  i <- 1
  if (fix_pi0) {
    opt_g$pi0 <- g$pi0
  } else {
    opt_g$pi0 <- 1 / (exp(-optpar[i]) + 1)
    i <- i + 1
  }

  if (fix_s2) {
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

  retlist <- opt_g
  retlist$val <- pn_llik_from_optval(optval, n1, n2, s2)

  return(retlist)
}

pn_llik_from_optval <- function(optval, n1, n2, s2) {
  if (length(s2) == 1) {
    sum.log.s2 <- n2 * log(s2)
  } else {
    sum.log.s2 <- sum(log(s2))
  }
  return(-optval - 0.5 * ((n1 + n2) * log(2 * pi) + sum.log.s2))
}

pn_partog <- function(g) {

}
