# Functions to compute MLE under point-normal prior.

mle_point_normal_logscale_grad <- function(x, s, g, control, fix_pi0, fix_mu) {
  # Optimization functions.
  fn <- function(par) {
    return(-loglik_point_normal(x,
                                s,
                                w = w_from_par(par, g, fix_pi0),
                                a = a_from_par(par, fix_pi0),
                                mu = mu_from_par(par, g, fix_pi0, fix_mu)))
  }
  gr <- function(par) {
    ret <- grad_negloglik_logscale_point_normal(x,
                                                s,
                                                w = w_from_par(par, g, fix_pi0),
                                                a = a_from_par(par, fix_pi0),
                                                mu = mu_from_par(par, g, fix_pi0, fix_mu))
    # Remove values that aren't being optimized over.
    return(ret[c(!fix_pi0, TRUE, !fix_mu)])
  }

  # Initial values.
  startpar <- numeric(0)
  if (!fix_pi0) {
    if (!is.null(g$pi0) && g$pi0 > 0 && g$pi0 < 1) {
      startpar <- c(startpar, log(1 / g$pi0 - 1))
    } else {
      startpar <- c(startpar, 0) # default for -logit(pi0)
    }
  }
  if (!is.null(g$a)) {
    startpar <- c(startpar, log(g$a))
  } else {
    startpar <- c(startpar, -log(mean(x^2))) # default for log(a)
  }
  if (!fix_mu) {
    if (!is.null(g$mu)) {
      startpar <- c(startpar, g$mu)
    } else {
      startpar <- c(startpar, mean(x)) # default for mu
    }
  }

  uu <- optimize_it(startpar, fn, gr, control,
                    mle_point_normal_hilo(x, s, fix_pi0, fix_mu),
                    x, s)

  # Pull pi0, a, and mu out of the optim results.
  retlist <- list()
  if (fix_pi0) {
    retlist$pi0 <- g$pi0
    retlist$a <- exp(uu$par[1])
    if (fix_mu) {
      retlist$mu <- g$mu
    } else {
      retlist$mu <- uu$par[2]
    }
  } else {
    retlist$pi0 <- 1 / (1 + exp(uu$par[1]))
    retlist$a <- exp(uu$par[2])
    if (fix_mu) {
      retlist$mu <- g$mu
    } else {
      retlist$mu <- uu$par[3]
    }
  }
  retlist$val <- uu$value
  return(retlist)
}


# Get w, a, and pi0 from par (the variables over which we're optimizing).
#   The parameters are -logit(pi0), log(a), and mu (regular scale), but we
#   need to be careful about the indexing because pi0 and mu can be fixed.
w_from_par <- function(par, g, fix_pi0) {
  if (fix_pi0)
    return(1 - g$pi0)
  else
    return(1 - (1/(1 + exp(par[1]))))
}

a_from_par <- function(par, fix_pi0) {
  if (fix_pi0)
    return(exp(par[1]))
  else
    return(exp(par[2]))
}

mu_from_par <- function(par, g, fix_pi0, fix_mu) {
  if (fix_pi0 && !fix_mu)
    return(par[2])
  else if (!fix_mu)
    return(par[3])
  else
    return(g$mu)
}


#' @importFrom stats optim
#'
optimize_it <- function(startpar, fn, gr, control, hilo, x, s) {
  uu <- try(optim(startpar, fn, gr, method = "BFGS",
                  control = control),
            silent=TRUE)

  # If optimization fails, try again with some limits; this should not
  # really happen but in preliminary testing sometimes we see optim
  # complain of infinite values, possibly because of extreme values of
  # the parameters?

  if (class(uu) == "try-error") {
    uu <- try(optim(startpar, fn, gr, method = "L-BFGS-B",
                    lower = hilo$lo, upper = hilo$hi, control = control))
  }

  # If optimization fails twice, save debug information and give up.

  if (class(uu) == "try-error") {
    saveRDS(list(startpar = startpar, x = x, s = s, control = control),
            "temp_debug.RDS")
    stop(paste("optim failed to converge; debug information saved to",
               "temp_debug.RDS"))
  }

  return(uu)
}


# Upper and lower bounds for optim in case the first attempt at optimization
#   fails.
mle_point_normal_hilo <- function(x, s, fix_pi0, fix_mu) {
  maxvar <- max(x^2)

  minvar <- (min(s) / 10)^2
  if (minvar < 1e-8) {
    minvar <- 1e-8
  }

  # Bounds for log(a):
  lo <- -log(maxvar)
  hi <- -log(minvar)

  if (!fix_pi0) {
    n <- length(x)
    lo <- c(log(1 / n), lo)
    hi <- c(log(n), hi)
  }

  if (!fix_mu) {
    lo <- c(lo, min(x) - 3 * max(s) - 3 * max(abs(x)))
    hi <- c(hi, max(x) + 3 * max(s) + 3 * max(abs(x)))
  }

  return(list(lo = lo, hi = hi))
}
