# @title Compute MLE under point-normal prior
#
# @description Computes MLE of w and a for data (x, s) under a
#   point-normal prior.
#
# @param x Observations.
#
# @param s Standard deviations.
#
# @param g List with initial values for g$a and g$pi0.
#
# @param control List of parameters to be passed to \code{optim}.
#
mle_point_normal_logscale_grad <- function(x, s, g, control) {
  # Do optimization of parameters on log scale (the parameters are
  #   -logit(pi0) and log(a)).
  
  maxvar <- max(x^2 - s^2) # Get upper bound on variance estimate.
  
  # Deal with case where everything is smaller than expected (null):
  if (maxvar < 0) {
    return(list(pi0 = 1, a = 1)) # Note that a is irrelevant if pi0 = 1.
  }
  
  # Set default starting point. This point is chosen based on the model
  #   where there is a single non-null value, based on the intuition that
  #   this is the case that is "hardest" to get right.
  #
  # if (is.null(startpar)) {
  #   startpar  <- c(log(1/length(x)),-log(maxvar))
  # }
  
  fn <- function(par) {
    -loglik_point_normal(x,
                         s,
                         w = 1 - (1/(1 + exp(par[1]))),
                         a = exp(par[2]))
  }
  
  gr <- function(par) {
    grad_negloglik_logscale_point_normal(x,
                                         s,
                                         w = 1 - (1/(1 + exp(par[1]))),
                                         a = exp(par[2]))
  }
  
  startpar <- c(0, 0) # default
  if (!is.null(g$pi0)) {
    startpar[1] <- log(1 / g$pi0 - 1)
  }
  if (!is.null(g$a)) {
    startpar[2] <- log(g$a)
  }
  
  uu <- optimize_it(startpar, fn, gr, control,
                    mle_point_normal_hilo(x, s, fix_pi0 = FALSE))
  
  return(list(pi0 = 1 / (1 + exp(uu$par[1])),
              a = exp(uu$par[2]),
              val = uu$value))
}


mle_point_normal_logscale_fixed_pi0 <- function(x, s, g, control) {
  
  fn <- function(par) {
    -loglik_point_normal(x, s, 1 - g$pi0, a = exp(par[1]))
  }
  
  gr <- function(par) {
    grad_negloglik_logscale_point_normal(x, s, 1 - g$pi0, a = exp(par[1]))[2]
  }
  
  if (!is.null(g$a)) {
    startpar <- log(g$a)
  } else {
    startpar <- 0 # default
  }
  
  uu <- optimize_it(startpar, fn, gr, control,
                    mle_point_normal_hilo(x, s, fix_pi0 = TRUE), x, s)
  
  return(list(pi0 = g$pi0, a = exp(uu$par[1]), val = uu$value))
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
  
  uu <- try(optim(startpar, fn, gr, method = "L-BFGS-B",
                  lower = hilo$lo, upper = hilo$hi, control = control))
  
  if (class(uu) == "try-error") {
    saveRDS(list(startpar = startpar, x = x, s = s, control = control),
            "temp_debug.RDS")
    stop(paste("optim failed to converge; debug information saved to",
               "temp_debug.RDS"))
  }
  
  return(uu)
}


# Get upper and lower bounds for optim in case the first attempt at
# optimization fails.
#
mle_point_normal_hilo <- function(x, s, fix_pi0) {
  maxvar <- max(x^2 - s^2)
  
  minvar <- (min(s) / 10)^2
  if (minvar < 1e-8) {
    minvar <- 1e-8
  }
  
  # Bounds for log(a):
  lo <- -log(maxvar)
  hi <- -log(minvar)
  
  if (!fix_pi0) {
    n <- length(x)
    lo <- c(log(1/n), lo)
    hi <- c(log(n), hi)
  }
  
  return(list(lo = lo, hi = hi))
}
