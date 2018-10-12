# @title Compute MLE under normal prior with no point-mass and arbitrary mean mu
#
# @description Computes MLE of mu and a for data (x, s) under a
#   normal prior.
#
# @param x Observations.
#
# @param s Standard deviations.
#
# @param g List with initial values for g$a and g$mu.
#
# @param control List of parameters to be passed to \code{optim}.
#
mle_normal_logscale_grad <- function(x, s, g, control) {
  # Do optimization of parameters on log scale (the parameters are
  #   mu and log(a)).
  
  # Set default starting point. This point is chosen based on the model
  #   where there is a single non-null value, based on the intuition that
  #   this is the case that is "hardest" to get right.
  #
  # if (is.null(startpar)) {
  #   startpar  <- c(log(1/length(x)),-log(maxvar))
  # }
  
  fn <- function(par) {
    -loglik_normal(x,
                   s,
                   mu = par[1],
                   a = exp(par[2]))
  }
  
  gr <- function(par) {
    grad_negloglik_logscale_normal(x,
                                   s,
                                   mu = par[1],
                                   a = exp(par[2]))
  }
  
  startpar <- c(mean(x), -log(max(1, var(x) - mean(s^2)))) # default
  if (!is.null(g$mu)) {
    startpar[1] <- g$mu
  }
  if (!is.null(g$a)) {
    startpar[2] <- log(g$a)
  }
  
  uu <- optimize_it_mu(startpar, fn, gr, control,
                    mle_point_normal_hilo(x, s, fix_mu = FALSE))
  
  return(list(mu = uu$par[1],
              a = exp(uu$par[2]),
              val = uu$value))
}


mle_normal_logscale_fixed_mu <- function(x, s, g, control) {
  
  fn <- function(par) {
    -loglik_normal(x, s, g$mu, a = exp(par[1]))
  }
  
  gr <- function(par) {
    grad_negloglik_logscale_normal(x, s, g$mu, a = exp(par[1]))[2]
  }
  
  if (!is.null(g$a)) {
    startpar <- log(g$a)
  } else {
    startpar <- -log(max(1, var(x) - mean(s^2))) # default
  }
  
  uu <- optimize_it(startpar, fn, gr, control,
                    mle_point_normal_hilo(x, s, fix_mu = TRUE), x, s)
  
  return(list(mu = g$mu, a = exp(uu$par[1]), val = uu$value))
}


#' @importFrom stats optim
#'
optimize_it_mu <- function(startpar, fn, gr, control, hilo, x, s) {
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
  
  if (class(uu) == "try-error") {
    saveRDS(list(startpar = startpar, x = x, s = s, control = control),
            "temp_debug.RDS")
    stop(paste("optim failed to converge; debug information saved to",
               "temp_debug.RDS"))
  }
  
  # check corner case where sigma^2 = 0 <=> log(a) = Inf, mu = sum(x_j / s_j^2) / sum(1 / s_j^2)
  mu_0 <- weighted.mean(x, 1/s^2) # solution for mu when sigma^2 = 0
  val_0 <- -loglik_normal(x, s, mu_0, Inf) # negloglik at this corner case
  
  if (val_0 < uu$value) { # if corner solution has smaller negloglik than optimizaed solution
    uu$par <- c(mu_0, Inf)
    uu$value <- val_0
  }
  
  return(uu)
}


# Get upper and lower bounds for optim in case the first attempt at
# optimization fails.
#
mle_normal_hilo <- function(x, s, fix_mu) {
  maxvar <- max(x^2 - s^2)
  
  minvar <- (min(s) / 10)^2
  if (minvar < 1e-8) {
    minvar <- 1e-8
  }
  
  # Bounds for log(a):
  lo <- -log(maxvar)
  hi <- -log(minvar)
  
  if (!fix_mu) {
    n <- length(x)
    lo <- c(min(x) - 6*max(s), lo)
    hi <- c(max(x) + 6*max(s), hi)
  }
  
  return(list(lo = lo, hi = hi))
}
