#' Constructor for horseshoe class
#'
#' Creates a horseshoe prior (see Carvalho, Polson, and Scott (2010)). The
#'   horseshoe is usually parametrized as
#'   \eqn{\theta_i \sim N(0, s^2 \tau^2 \lambda_i^2)},
#'   \eqn{\lambda_i \sim \mathrm{Cauchy}^+(0, 1)},
#'   with \eqn{s^2} the variance of the error distribution. We use a single
#'   parameter \code{scale}, which corresponds to \eqn{s\tau} and thus does
#'   not depend on the error distribution.
#'
#' @param scale The scale parameter (must be a scalar).
#'
#' @return An object of class \code{horseshoe} (a list with a single element
#'   \code{scale}, described above).
#'
#' @export
#'
horseshoe <- function(scale) {
  structure(data.frame(scale), class = "horseshoe")
}

#' @importFrom horseshoe HS.post.mean HS.post.var HS.normal.means
#'
horseshoe_workhorse <- function(x = x,
                                s = s,
                                mode = mode,
                                scale = scale,
                                g_init = g_init,
                                fix_g = fix_g,
                                output = output,
                                control = control,
                                call = call) {
  if (length(s) != 1) {
    if (!isTRUE(all.equal(min(s) / mean(s), max(s) / mean(s)))) {
      stop("The horseshoe prior family requires homoskedastic standard errors.")
    } else {
      s <- s[1]
    }
  }

  if (mode != 0) {
    stop("Nonzero modes have not yet been implemented for the horseshoe prior family.")
  }
  call$mode <- NULL # Bypasses warning about setting mode = 0 when g is fixed.

  check_g_init(g_init = g_init,
               fix_g = fix_g,
               mode = mode,
               scale = scale,
               pointmass = FALSE,
               call = call,
               class_name = "horseshoe",
               scale_name = "scale")

  if (identical(scale, "estimate")) {
    optres <- myHS.MMLE(y = x, Sigma2 = s^2, control = control)
    tau <- optres$maximum
    obj <- optres$objective
  } else {
    tau <- scale / s
    obj <- MMLE.M(tau, data = x, data.var = s^2)
  }

  retlist <- list()

  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, x, s)
  }

  if (posterior_in_output(output)) {
    posterior <- list()

    if (result_in_output(output)) {
      posterior$mean  <- horseshoe::HS.post.mean(x, tau, s^2)
      posterior$sd    <- sqrt(horseshoe::HS.post.var(x, tau, s^2))
      posterior$mean2 <- posterior$mean^2 + posterior$sd^2
    }

    if (lfsr_in_output(output)) {
      warning("LFSR calculations are not implemented for the horseshoe prior",
              " family. Please use the posterior sampler instead.")
    }

    retlist <- add_posterior_to_retlist(retlist, posterior, output)
  }

  if (g_in_output(output)) {
    fitted_g <- horseshoe(tau * s)
    retlist  <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    n       <- length(x)
    loglik  <- obj - 1.5 * n * log(pi) - 0.5 * n * log(2) - n * sum(log(s))
    retlist <- add_llik_to_retlist(retlist, loglik)
  }

  if (sampler_in_output(output)) {
    post_sampler <- function(nsamp, burn = 1000) {
      cat("MCMC Sampling with", burn, "burn-in samples\n")
      samp <- HS.normal.means(x, tau = tau, Sigma2 = s^2, nmc = nsamp, burn = burn)
      return(samp$BetaSamples)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(retlist)
}


# A lightly modified version of HS.MMLE from package horseshoe (version 0.2.0):

#' @importFrom stats optimize integrate
#'
myHS.MMLE <- function(y, Sigma2, control){
  do.call(optimize, c(list(f = MMLE.M,
                           data = y,
                           data.var = Sigma2,
                           lower = (1 / length(y)),
                           upper = 1,
                           maximum = TRUE),
                      control))
}


# Following functions are copy-pasted from package horseshoe (version 0.2.0):

## Helper functions for estimating the MMLE
## Can handle very sparse situations only for n <= 425
Basic.MMLE <- function(u, y, t, data.var){ #integrand
  (1-u)^(-1/2)*(1 - (1-t^2)*u)^(-1)*exp(-u*y^2 / (2*data.var))
}

## Helper functions for estimating the MMLE
## Can handle very sparse situations only for n <= 425
Term.MMLE <- function(y, t, data.var){ #integrate the integrand and take the log
  log( stats::integrate(Basic.MMLE, 0, 1, y = y,  t = t, data.var = data.var, rel.tol = 1e-8)$value ) + log(t)
}

## Helper functions for estimating the MMLE
## Can handle very sparse situations only for n <= 425
Term.MMLE.vec <- Vectorize(Term.MMLE) #handles a vector as input

## Helper functions for estimating the MMLE
## Can handle very sparse situations only for n <= 425
#' @keywords internal
MMLE.M <- function(t, data, data.var){ #the quantity to be maximized
  sum(Term.MMLE.vec(data, t, data.var))
}
