# The function that nlm optimizes when fitting a point-Laplace distribution.

#' @importFrom stats pnorm
#'
pl_nllik <- function(par,
                     fix_pi0, fix_a, fix_mu, alpha, beta, mu,
                     x, s,
                     calc_grad, calc_hess) {
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
  }

  w <- 1 / (1 + exp(-alpha))
  a <- exp(beta)

  # Write the negative log likelihood as -log((1 - w)f + wg), where f
  #   corresponds to the point mass and g to the Laplace component.

  # Point mass:
  lf <- -0.5 * log(2 * pi * s^2) - 0.5 * (x - mu)^2 / s^2

  # This part comes from the left tail.
  xleft <- (x - mu) / s + s * a
  lpnormleft <- pnorm(xleft, log.p = TRUE, lower.tail = FALSE)
  lgleft <- log(a / 2) + s^2 * a^2 / 2 + a * (x - mu) + lpnormleft

  # This part comes from the right tail.
  xright <- (x - mu) / s - s * a
  lpnormright <- pnorm(xright, log.p = TRUE)
  lgright <- log(a / 2) + s^2 * a^2 / 2 - a * (x - mu) + lpnormright

  # Laplace component:
  lg <- logscale_add(lgleft, lgright)

  llik <- logscale_add(log(1 - w) + lf, log(w) + lg)
  nllik <- -sum(llik)

  # Gradients:
  if (calc_grad || calc_hess) {
    # Since the likelihood appears in the denominator of all gradients, all
    #   quantities below divide by the likelihood (but this is not reflected in
    #   the variable names!).
    grad <- numeric(length(par))
    i <- 1
    if (!fix_pi0) {
      # Derivatives with respect to w and alpha (= par[1]).
      f <- exp(lf - llik)
      g <- exp(lg - llik)
      dnllik.dw <- f - g
      dw.dalpha <- w * (1 - w)
      dnllik.dalpha <- dnllik.dw * dw.dalpha

      grad[i] <- sum(dnllik.dalpha)
      i <- i + 1
    }
    if (!fix_a || !fix_mu) {
      dlogpnorm.left <- -exp(-log(2 * pi) / 2 - xleft^2 / 2 - lpnormleft)
      dlogpnorm.right <- exp(-log(2 * pi) / 2 - xright^2 / 2 - lpnormright)
    }
    if (!fix_a) {
      # Derivatives with respect to a and beta (= par[2]).
      dgleft.da <- exp(lgleft - llik) * (1 / a + a * s^2 + (x - mu) + s * dlogpnorm.left)
      dgright.da <- exp(lgright - llik) * (1 / a + a * s^2 - (x - mu) - s * dlogpnorm.right)
      dg.da <- dgleft.da + dgright.da
      dnllik.da <- -w * dg.da
      da.dbeta <- a
      dnllik.dbeta <- dnllik.da * da.dbeta

      grad[i] <- sum(dnllik.dbeta)
      i <- i + 1
    }
    if (!fix_mu) {
      # Derivatives with respect to mu.
      df.dmu <- exp(lf - llik) * ((x - mu) / s^2)
      dgleft.dmu <- exp(lgleft - llik) * (-a - dlogpnorm.left / s)
      dgright.dmu <- exp(lgright - llik) * (a - dlogpnorm.right / s)
      dg.dmu <- dgleft.dmu + dgright.dmu
      dnllik.dmu <- -(1 - w) * df.dmu - w * dg.dmu

      grad[i] <- sum(dnllik.dmu)
    }

    attr(nllik, "gradient") <- grad
  }

  if (calc_hess) {
    hess <- matrix(nrow = length(par), ncol = length(par))
    i <- 1
    if (!fix_pi0) {
      # Derivatives with respect to w and alpha (= par[1]).
      # Second derivative with respect to w (alpha).
      d2nllik.dw2 <- (dnllik.dw)^2
      d2w.dalpha2 <- (1 - 2 * w) * dw.dalpha
      d2nllik.dalpha2 <- d2nllik.dw2 * (dw.dalpha)^2 + dnllik.dw * (d2w.dalpha2)

      hess[i, i] <- sum(d2nllik.dalpha2)

      j <- i + 1
      if (!fix_a) {
        # Mixed derivative with respect to w and a (alpha and beta).
        d2nllik.dwda <- dnllik.dw * dnllik.da - dg.da
        d2nllik.dalphadbeta <- d2nllik.dwda * dw.dalpha * da.dbeta

        hess[i, j] <- hess[j, i] <- sum(d2nllik.dalphadbeta)
        j <- j + 1
      }
      if (!fix_mu) {
        # Mixed derivative with respect to w (alpha) and mu.
        d2nllik.dwdmu <- dnllik.dw * dnllik.dmu - (dg.dmu - df.dmu)
        d2nllik.dalphadmu <- d2nllik.dwdmu * dw.dalpha

        hess[i, j] <- hess[j, i] <- sum(d2nllik.dalphadmu)
      }

      i <- i + 1
    }
    if (!fix_a) {
      # Second derivative with respect to a (beta).
      d2gleft.da2 <- dgleft.da * (1 / a + a * s^2 + (x - mu) + s * dlogpnorm.left) +
        exp(lgleft - llik) * (-1 / a^2 + s^2 * (1 - dlogpnorm.left * xleft - dlogpnorm.left^2))
      d2gright.da2 <- dgright.da * (1 / a + a * s^2 - (x - mu) - s * dlogpnorm.right) +
        exp(lgright - llik) * (-1 / a^2 + s^2 * (1 - dlogpnorm.right * xright - dlogpnorm.right^2))
      d2g.da2 <- d2gleft.da2 + d2gright.da2
      d2nllik.da2 <- (dnllik.da)^2 - w * d2g.da2
      d2a.dbeta2 <- da.dbeta
      d2nllik.dbeta2 <- d2nllik.da2 * (da.dbeta)^2 + dnllik.da * (d2a.dbeta2)

      hess[i, i] <- sum(d2nllik.dbeta2)

      j <- i + 1
      if (!fix_mu) {
        # Mixed derivative with respect to a (beta) and mu.
        d2gleft.dadmu <- dgleft.da * (-a - dlogpnorm.left / s) +
          exp(lgleft - llik) * (-1 + dlogpnorm.left * xleft + dlogpnorm.left^2)
        d2gright.dadmu <- dgright.da * (a - dlogpnorm.right / s) +
          exp(lgright - llik) * (1 - dlogpnorm.right * xright - dlogpnorm.right^2)
        d2g.dadmu <- d2gleft.dadmu + d2gright.dadmu
        d2nllik.dbetadmu <- dnllik.dbeta * dnllik.dmu - w * d2g.dadmu * da.dbeta

        hess[i, j] <- hess[j, i] <- sum(d2nllik.dbetadmu)
      }

      i <- i + 1
    }
    if (!fix_mu) {
      # Second derivative with respect to mu.
      d2f.dmu2 <- df.dmu * ((x - mu) / s^2) - exp(lf - llik) / s^2
      d2gleft.dmu2 <- dgleft.dmu * (-a - dlogpnorm.left / s) -
        exp(lgleft - llik) * (dlogpnorm.left * xleft + dlogpnorm.left^2) / s^2
      d2gright.dmu2 <- dgright.dmu * (a - dlogpnorm.right / s) -
        exp(lgright - llik) * (dlogpnorm.right * xright + dlogpnorm.right^2) / s^2
      d2g.dmu2 <- d2gleft.dmu2 + d2gright.dmu2
      d2nllik.dmu2 <- (dnllik.dmu)^2 - (1 - w) * d2f.dmu2 - w * d2g.dmu2

      hess[i, i] <- sum(d2nllik.dmu2)
    }

    attr(nllik, "hessian") <- hess
  }

  return(nllik)
}

logscale_add <- function(log.x, log.y) {
  C <- pmax(log.x, log.y)
  return(log(exp(log.x - C) + exp(log.y - C)) + C)
}
