#' Constructor for gammamix class
#'
#' Creates a finite mixture of gamma distributions.
#'
#' @param pi A vector of mixture proportions.
#'
#' @param shape A vector of shape parameters.
#'
#' @param scale A vector of scale parameters.
#'
#' @param shift A vector of shift parameters.
#'
#' @return An object of class \code{gammamix} (a list with elements
#'   \code{pi}, \code{shape}, \code{scale}, and \code{shift}, described above).
#'
#' @export
#'
gammamix <- function(pi, shape, scale, shift = rep(0, length(pi))) {
  structure(data.frame(pi, shape, scale, shift), class = "gammamix")
}

#' @importFrom stats pgamma
#' @importFrom ashr comp_cdf
#'
#' @method comp_cdf gammamix
#'
#' @export
#'
comp_cdf.gammamix = function (m, y, lower.tail = TRUE) {
  pshiftgamma <- function(q, shape, scale = 1, shift = 0, lower.tail = TRUE) {
    q <- q - shift
    scale <- ifelse(scale == 0, 1e-16, scale)
    return(pgamma(q, shape, scale = scale, lower.tail = lower.tail))
  }
  return(vapply(y, pshiftgamma, m$shift, m$shape, m$scale, m$shift, lower.tail))
}

# The point-exponential family uses the above ebnm class gammamix.
#
pe_checkg <- function(g_init, fix_g, mode, scale, pointmass, call) {
  check_g_init(g_init = g_init,
               fix_g = fix_g,
               mode = mode,
               scale = scale,
               pointmass = pointmass,
               call = call,
               class_name = "gammamix",
               scale_name = "scale",
               mode_name = "shift")

  if (!is.null(g_init)
      && !isTRUE(all.equal(g_init$shape, rep(1, length(g_init$shape))))) {
    stop("g_init must be of class gammamix with shape parameter = 1")
  }
}


# Point-exponential parameters are alpha = -logit(pi0), beta = log(a), and mu.
#
pe_initpar <- function(g_init, mode, scale, pointmass, x, s) {
  if (!is.null(g_init) && length(g_init$pi) == 1) {
    par <- list(alpha = Inf,
                beta = -log(g_init$scale),
                mu = g_init$shift)
  } else if (!is.null(g_init) && length(g_init$pi) == 2) {
    par <- list(alpha = log(1 / g_init$pi[1] - 1),
                beta = -log(g_init$scale[2]),
                mu = g_init$shift[1])
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
      par$mu <- min(x) # default
    }
  }

  return(par)
}


pe_scalepar <- function(par, scale_factor) {
  if (!is.null(par$beta)) {
    par$beta <- par$beta - log(scale_factor)
  }
  if (!is.null(par$mu)) {
    par$mu <- scale_factor * par$mu
  }

  return(par)
}


# No precomputations are done for point-exponential.
#
pe_precomp <- function(x, s, par_init, fix_par) {
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
pe_nllik <- function(par, x, s, par_init, fix_par,
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
  #   corresponds to the point mass and g to the exponential component.

  # Point mass:
  lf <- -0.5 * log(2 * pi * s^2) - 0.5 * (x - mu)^2 / s^2

  # Exponential component:
  xright <- (x - mu) / s - s * a
  lpnormright <- pnorm(xright, log.p = TRUE)
  lg <- log(a) + s^2 * a^2 / 2 - a * (x - mu) + lpnormright

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
      dlogpnorm.right <- exp(-log(2 * pi) / 2 - xright^2 / 2 - lpnormright)
    }
    if (!fix_a) {
      # Derivatives with respect to a and beta (= par[2]).
      dg.da <- exp(lg - llik) * (1 / a + a * s^2 - (x - mu) - s * dlogpnorm.right)
      dnllik.da <- -w * dg.da
      da.dbeta <- a
      dnllik.dbeta <- dnllik.da * da.dbeta

      grad[i] <- sum(dnllik.dbeta)
      i <- i + 1
    }
    if (!fix_mu) {
      # Derivatives with respect to mu.
      df.dmu <- exp(lf - llik) * ((x - mu) / s^2)
      dg.dmu <- exp(lg - llik) * (a - dlogpnorm.right / s)
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
      d2g.da2 <- dg.da * (1 / a + a * s^2 - (x - mu) - s * dlogpnorm.right) +
        exp(lg - llik) * (-1 / a^2 + s^2 * (1 - dlogpnorm.right * xright - dlogpnorm.right^2))
      d2nllik.da2 <- (dnllik.da)^2 - w * d2g.da2
      d2a.dbeta2 <- da.dbeta
      d2nllik.dbeta2 <- d2nllik.da2 * (da.dbeta)^2 + dnllik.da * (d2a.dbeta2)

      hess[i, i] <- sum(d2nllik.dbeta2)

      j <- i + 1
      if (!fix_mu) {
        # Mixed derivative with respect to a (beta) and mu.
        d2g.dadmu <- dg.da * (a - dlogpnorm.right / s) +
          exp(lg - llik) * (1 - dlogpnorm.right * xright - dlogpnorm.right^2)
        d2nllik.dbetadmu <- dnllik.dbeta * dnllik.dmu - w * d2g.dadmu * da.dbeta

        hess[i, j] <- hess[j, i] <- sum(d2nllik.dbetadmu)
      }

      i <- i + 1
    }
    if (!fix_mu) {
      # Second derivative with respect to mu.
      d2f.dmu2 <- df.dmu * ((x - mu) / s^2) - exp(lf - llik) / s^2
      d2g.dmu2 <- dg.dmu * (a - dlogpnorm.right / s) -
        exp(lg - llik) * (dlogpnorm.right * xright + dlogpnorm.right^2) / s^2
      d2nllik.dmu2 <- (dnllik.dmu)^2 - (1 - w) * d2f.dmu2 - w * d2g.dmu2

      hess[i, i] <- sum(d2nllik.dmu2)
    }

    attr(nllik, "hessian") <- hess
  }

  return(nllik)
}


# Postcomputations: check boundary solutions.
#
pe_postcomp <- function(optpar, optval, x, s, par_init, fix_par, scale_factor) {
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
pe_summres <- function(x, s, optpar, output) {
  w  <- 1 - 1 / (exp(optpar$alpha) + 1)
  a  <- exp(optpar$beta)
  mu <- optpar$mu

  return(pe_summres_untransformed(x, s, w, a, mu, output))
}

#' @importFrom ashr my_etruncnorm my_e2truncnorm
#'
pe_summres_untransformed <- function(x, s, w, a, mu, output) {
  x <- x - mu

  wpost <- wpost_exp(x, s, w, a)

  post <- list()

  if (result_in_output(output)) {
    post$mean  <- wpost * my_etruncnorm(0, Inf, x - s^2 * a, s)
    post$mean2 <- wpost * my_e2truncnorm(0, Inf, x - s^2 * a, s)

    if (any(is.infinite(s))) {
      post$mean[is.infinite(s)]  <- w / a
      post$mean2[is.infinite(s)] <- 2 * w / a^2
    }
    post$sd <- sqrt(pmax(0, post$mean2 - post$mean^2))

    post$mean2 <- post$mean2 + mu^2 + 2 * mu * post$mean
    post$mean  <- post$mean + mu
  }

  if ("lfsr" %in% output) {
    post$lfsr <- 1 - wpost
    if (any(is.infinite(s))) {
      post$lfsr[is.infinite(s)] <- 1 - w
    }
  }

  return(post)
}

#  Calculate posterior weights for non-null effects.
#
#' @importFrom stats dnorm pnorm
#'
wpost_exp <- function(x, s, w, a) {
  if (w == 0) {
    return(rep(0, length(x)))
  }

  if (w == 1) {
    return(rep(1, length(x)))
  }

  lf <- dnorm(x, 0, s, log = TRUE)
  lg <- log(a) + s^2 * a^2 / 2 - a * x + pnorm(x / s - s * a, log.p = TRUE)
  wpost <- w / (w + (1 - w) * exp(lf - lg))

  return(wpost)
}


# Point-exponential parameters are alpha = -logit(pi0), beta = log(a), and mu.
#   The above ebnm class gammamix is used.
#
pe_partog <- function(par) {
  pi0   <- 1 / (exp(par$alpha) + 1)
  scale <- exp(-par$beta)
  mode  <- par$mu

  if (pi0 == 0) {
    g <- gammamix(pi = 1,
                  shape = 1,
                  scale = scale,
                  shift = mode)
  } else {
    g <- gammamix(pi = c(pi0, 1 - pi0),
                  shape = c(1, 1),
                  scale = c(0, scale),
                  shift = rep(mode, 2))
  }

  return(g)
}


# Sample from the posterior under point-exponential prior.
#
# @param nsamp The number of samples to return per observation.
#
# @return An nsamp by length(x) matrix containing samples from the
#   posterior, with each row corresponding to a single sample.
#
pe_postsamp <- function(x, s, optpar, nsamp) {
  w  <- 1 - 1 / (exp(optpar$alpha) + 1)
  a  <- exp(optpar$beta)
  mu <- optpar$mu

  return(pe_postsamp_untransformed(x, s, w, a, mu, nsamp))
}

#' @importFrom truncnorm rtruncnorm
#' @importFrom stats rbinom
#'
pe_postsamp_untransformed <- function(x, s, w, a, mu, nsamp) {
  x <- x - mu

  wpost <- wpost_exp(x, s, w, a)

  nobs <- length(wpost)

  is_nonnull <- rbinom(nsamp * nobs, 1, rep(wpost, each = nsamp))

  if (length(s) == 1) {
    s <- rep(s, nobs)
  }

  # rtruncnorm is not vectorized.
  positive_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, 0, Inf, mean, sd)
  }, mean = x - s^2 * a, sd = s)

  samp <- matrix(0, nrow = nsamp, ncol = length(wpost))
  samp[is_nonnull == 1] <- positive_samp[is_nonnull == 1]

  samp <- samp + mu

  return(samp)
}
