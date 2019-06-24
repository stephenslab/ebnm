# Computes MLE for g under point-normal prior.
#
#' @importFrom stats optim
#'
mle_point_laplace <- function(x, s, startpar = NULL, control = NULL) {
  startpar <- pl_startpar(x, s)

  lf <- -0.5 * log(2 * pi * s^2) - 0.5 * x^2 / s^2

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

  return(retlist)
}

# Initial values.
pl_startpar <- function(x, s) {
  # Defaults for -logit(pi0) and log(a).
  return(c(0, -log(mean(x^2))))
}

pl_g_from_optpar <- function(par) {
  return(list(pi0 = 1 / (1 + exp(par[1])), a = exp(par[2])))
}
