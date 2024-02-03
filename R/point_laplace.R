#' Constructor for laplacemix class
#'
#' Creates a finite mixture of Laplace distributions.
#'
#' @param pi A vector of mixture proportions.
#'
#' @param mean A vector of means.
#'
#' @param scale A vector of scale parameters.
#'
#' @return An object of class \code{laplacemix} (a list with elements
#'   \code{pi}, \code{mean}, and \code{scale}, described above).
#'
#' @export
#'
laplacemix <- function(pi, mean, scale) {
  structure(data.frame(pi, mean, scale), class="laplacemix")
}

#' @importFrom stats pexp
#' @importFrom ashr comp_cdf
#'
#' @method comp_cdf laplacemix
#'
#' @export
#'
comp_cdf.laplacemix = function (m, y, lower.tail = TRUE) {
  plaplace <- function(q, location = 0, rate = 1, lower.tail = TRUE) {
    q <- q - location
    if (!lower.tail) {
      q <- -q
    }
    return(ifelse(q < 0,
                  pexp(abs(q), rate, lower.tail = FALSE) / 2,
                  0.5 + pexp(abs(q), rate) / 2
    ))
  }
  return(vapply(y, plaplace, m$mean, m$mean, 1 / m$scale, lower.tail))
}

# The point-Laplace family uses the above ebnm class laplacemix.
#
pl_checkg <- function(g_init, fix_g, mode, scale, pointmass, call) {
  check_g_init(g_init = g_init,
               fix_g = fix_g,
               mode = mode,
               scale = scale,
               pointmass = pointmass,
               call = call,
               class_name = "laplacemix",
               scale_name = "scale")
}


# Point-Laplace parameters are alpha = -logit(pi0), beta = -log(lambda), and mu.
#
pl_initpar <- function(g_init, mode, scale, pointmass, x, s) {
  if (!is.null(g_init) && length(g_init$pi) == 1) {
    par <- list(alpha = Inf,
                beta = -log(g_init$scale),
                mu = g_init$mean)
  } else if (!is.null(g_init) && length(g_init$pi) == 2) {
    par <- list(alpha = log(1 / g_init$pi[1] - 1),
                beta = -log(g_init$scale[2]),
                mu = g_init$mean[1])
  } else {
    par <- list()
    if (!pointmass) {
      par$alpha <- Inf
    } else {
      par$alpha <- 0 # default
    }
    if (!identical(scale, "estimate")) {
      if (length(scale) != 1) {
        stop("Argument 'scale' must be either 'estimate' or a scalar.")
      }
      par$beta <- -log(scale)
    } else {
      par$beta <- -0.5 * log(mean(x^2) / 2) # default
    }
    if (!identical(mode, "estimate")) {
      par$mu <- mode
    } else {
      par$mu <- mean(x) # default
    }
  }

  return(par)
}


pl_scalepar <- function(par, scale_factor) {
  if (!is.null(par$beta)) {
    par$beta <- par$beta - log(scale_factor)
  }
  if (!is.null(par$mu)) {
    par$mu <- scale_factor * par$mu
  }

  return(par)
}


# No precomputations are done for point-Laplace.
#
pl_precomp <- function(x, s, par_init, fix_par) {
  fix_mu  <- fix_par[3]

  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }

  return(NULL)
}


# The negative log likelihood.
#
#' @importFrom stats pnorm
#'
pl_nllik <- function(par, x, s, par_init, fix_par,
                     calc_grad, calc_hess) {
  fix_pi0 <- fix_par[1]
  fix_a   <- fix_par[2]
  fix_mu  <- fix_par[3]

  p <- unlist(par_init)
  p[!fix_par] <- par

  w <- 1 - 1 / (1 + exp(p[1]))
  a <- exp(p[2])
  mu <- p[3]

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


# Postcomputations: check boundary solutions.
#
pl_postcomp <- function(optpar, optval, x, s, par_init, fix_par, scale_factor) {
  llik <- -optval
  retlist <- list(par = optpar, val = llik)

  # Check the solution pi0 = 1.
  fix_pi0 <- fix_par[1]
  fix_mu  <- fix_par[3]
  if (!fix_pi0 && fix_mu) {
    pi0_llik <- sum(-0.5 * log(2 * pi * s^2) - 0.5 * (x - par_init$mu)^2 / s^2)
    pi0_llik <- pi0_llik + sum(is.finite(x)) * log(scale_factor)
    if (pi0_llik > llik) {
      retlist$par$alpha <- -Inf
      retlist$par$beta <- 0
      retlist$val <- pi0_llik
    }
  }

  return(retlist)
}


# Summary results.
#
pl_summres <- function(x, s, optpar, output) {
  w  <- 1 - 1 / (exp(optpar$alpha) + 1)
  a  <- exp(optpar$beta)
  mu <- optpar$mu

  return(pl_summres_untransformed(x, s, w, a, mu, output))
}

#' @importFrom ashr my_etruncnorm my_e2truncnorm
#'
pl_summres_untransformed <- function(x, s, w, a, mu, output) {
  x <- x - mu

  wpost <- wpost_laplace(x, s, w, a)
  lm <- lambda(x, s, a)

  post <- list()

  if (result_in_output(output)) {
    post$mean  <- wpost * (lm * my_etruncnorm(0, Inf, x - s^2 * a, s)
                           + (1 - lm) * my_etruncnorm(-Inf, 0, x + s^2 * a, s))
    post$mean2 <- wpost * (lm * my_e2truncnorm(0, Inf, x - s^2 * a, s)
                           + (1 - lm) * my_e2truncnorm(-Inf, 0, x + s^2 * a, s))
    if (any(is.infinite(s))) {
      post$mean[is.infinite(s)]  <- 0
      post$mean2[is.infinite(s)] <- 2 * w / a^2
    }
    post$mean2 <- pmax(post$mean2, post$mean^2)
    post$sd    <- sqrt(pmax(0, post$mean2 - post$mean^2))

    post$mean2 <- post$mean2 + mu^2 + 2 * mu * post$mean
    post$mean  <- post$mean + mu
  }

  if ("lfsr" %in% output) {
    post$lfsr <- (1 - wpost) + wpost * pmin(lm, 1 - lm)
    if (any(is.infinite(s))) {
      post$lfsr[is.infinite(s)] <- 1 - w / 2
    }
  }

  return(post)
}

#  Calculate posterior weights for non-null effects.
#
#' @importFrom stats dnorm
#'
wpost_laplace <- function(x, s, w, a) {
  if (w == 0) {
    return(rep(0, length(x)))
  }

  if (w == 1) {
    return(rep(1, length(x)))
  }

  lf <- dnorm(x, 0, s, log = TRUE)
  lg <- logg_laplace(x, s, a)
  wpost <- w / (w + (1 - w) * exp(lf - lg))

  return(wpost)
}

# This is the log of g, Laplace(a) convolved with normal.
#
#' @importFrom stats pnorm
#'
logg_laplace = function(x, s, a) {
  lg1 <- -a * x + pnorm((x - s^2 * a) / s, log.p = TRUE)
  lg2 <-  a * x + pnorm((x + s^2 * a) / s, log.p = TRUE, lower.tail = FALSE)
  lfac <- pmax(lg1, lg2)
  return(log(a / 2) + s^2 * a^2 / 2 + lfac + log(exp(lg1 - lfac) + exp(lg2 - lfac)))
}

# Compute the lambda function equation (2.7) from Kan Xu's thesis, which is
#   the posterior probability of being positive given a non-zero effect.
#
#' @importFrom stats pnorm
#'
lambda <- function(x, s, a) {
  lm1 <- -a * x + pnorm(x / s - s * a, log.p = TRUE)
  lm2 <-  a * x + pnorm(x / s + s * a, log.p = TRUE, lower.tail = FALSE)

  lm <- 1 / (1 + exp(lm2 - lm1))

  return(lm)
}


# Point-Laplace parameters are alpha = -logit(pi0), beta = -log(lambda), and mu.
#   The above ebnm class laplacemix is used.
#
pl_partog <- function(par) {
  pi0   <- 1 / (exp(par$alpha) + 1)
  scale <- exp(-par$beta)
  mean  <- par$mu

  if (pi0 == 0) {
    g <- laplacemix(pi = 1,
                    mean = mean,
                    scale = scale)
  } else {
    g <- laplacemix(pi = c(pi0, 1 - pi0),
                    mean = rep(mean, 2),
                    scale = c(0, scale))
  }

  return(g)
}


# Sample from the posterior under point-Laplace prior.
#
# @param nsamp The number of samples to return per observation.
#
# @return An nsamp by length(x) matrix containing samples from the
#   posterior, with each row corresponding to a single sample.
#
pl_postsamp <- function(x, s, optpar, nsamp) {
  w  <- 1 - 1 / (exp(optpar$alpha) + 1)
  a  <- exp(optpar$beta)
  mu <- optpar$mu

  return(pl_postsamp_untransformed(x, s, w, a, mu, nsamp))
}

#' @importFrom truncnorm rtruncnorm
#' @importFrom stats rbinom
#'
pl_postsamp_untransformed <- function(x, s, w, a, mu, nsamp) {
  x <- x - mu

  wpost <- wpost_laplace(x, s, w, a)
  l <- lambda(x, s, a)

  nobs <- length(wpost)

  is_nonnull <- rbinom(nsamp * nobs, 1, rep(wpost, each = nsamp))
  is_positive <- rbinom(nsamp * nobs, 1, rep(l, each = nsamp))

  if (length(s) == 1) {
    s <- rep(s, nobs)
  }

  # rtruncnorm is not vectorized.
  negative_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, -Inf, 0, mean, sd)
  }, mean = x + s^2 * a, sd = s)
  positive_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, 0, Inf, mean, sd)
  }, mean = x - s^2 * a, sd = s)

  samp <- matrix(0, nrow = nsamp, ncol = length(wpost))
  samp[is_nonnull & is_positive] <- positive_samp[is_nonnull & is_positive]
  samp[is_nonnull & !is_positive] <- negative_samp[is_nonnull & !is_positive]

  samp <- samp + mu

  return(samp)
}
