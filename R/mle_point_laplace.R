# Computes MLE for g under point-normal prior.
#
#' @importFrom stats optim
#'
mle_point_laplace <- function(x, s, g, control, fix_a) {
  startpar <- pl_startpar(x, s, g)

  lf <- calc_lf(x, s)

  fn_params <- list(x = x, s = s, lf = lf)

  optres <- try(do.call(nlm, c(list(pl_nlm_fn, startpar), fn_params, control)),
                silent = TRUE)

  if (inherits(optres, "try-error")) {
    warning("First optimization attempt failed. Retrying with fewer ",
            "significant digits.")
    control <- modifyList(control, list(ndigit = 8))
    optres <- try(do.call(nlm, c(list(pl_nlm_fn, startpar), fn_params, control)),
                  silent = TRUE)
  }

  if (inherits(optres, "try-error")) {
    stop("Re-attempt at optimization failed, possibly due to one or more ",
         "very small standard errors. The smallest nonzero standard error ",
         "seen during optimization was ", signif(min(s), 2), ".")
  }

  retlist <- pl_g_from_optpar(optres$estimate)
  retlist$val <- -optres$minimum

  # Check the solution pi0 = 1.
  if (sum(lf) > retlist$val) {
    retlist$pi0 <- 1
    retlist$a <- 1
    retlist$val <- sum(lf)
  }

  return(retlist)
}

# Initial values.
pl_startpar <- function(x, s, g) {
  startpar <- numeric(0)

  if (!is.null(g$pi0) && g$pi0 > 0 && g$pi0 < 1) {
    startpar <- c(startpar, log(1 / g$pi0 - 1))
  } else {
    startpar <- c(startpar, 0) # default for -logit(pi0)
  }

  if (!is.null(g$a)) {
    startpar <- c(startpar, log(g$a))
  } else {
    startpar <- c(startpar, -0.5 * log(mean(x^2) / 2)) # default for log(a)
  }

  return(startpar)
}

pl_g_from_optpar <- function(par) {
  return(list(pi0 = 1 / (1 + exp(par[1])), a = exp(par[2])))
}

calc_lf <- function(x, s) {
  return(-0.5 * log(2 * pi * s^2) - 0.5 * x^2 / s^2)
}
