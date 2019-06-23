#' Solve the EBNM problem
#'
#' Solves the empirical Bayes normal means problem using a specified class of
#'   priors.
#'
#' @details TODO: update me.
#'
#' Given vectors of data \code{x} and standard errors \code{s},
#'   solve the EBNM problem with a point-normal or point-laplace prior. The
#'   model is \deqn{x_j \sim N(\theta_j, s_j^2),} where \eqn{s_j} are given and
#'   \eqn{\theta_j \sim g}, with \eqn{g} either a mixture of a point mass at
#'   \eqn{\mu} and a normal distribution: \deqn{\theta_j \sim \pi_0 \delta_\mu
#'   + (1 - \pi_0)N(\mu, 1/a)} or a mixture of a point mass at zero and a
#'   laplace distribution: \deqn{\theta_j \sim \pi_0 \delta_0 +
#'   (1 - \pi_0)DExp(a).} \eqn{\pi_0}, \eqn{a}, and \eqn{\mu} are estimated by
#'   marginal maximum likelihood.
#'
#' @param x A vector of observations.
#'
#' @param s A vector of standard deviations (or a scalar if all are equal).
#'
#' @param mode Fixes the location of the prior mode. Set to "estimate" to
#'   estimate it from the data.
#'
#' @param scale Fixes the scale of the prior. Corresponds to the standard
#'   deviation of the normal component for normal and point-normal
#'   distributions; the rate parameter of the Laplace component for
#'   point-Laplace distributions; and parameter \code{mixsd} for adaptive
#'   shrinkage priors (see \code{\link[ashr]{ash}}). Set to "estimate" to
#'   estimate it from the data.
#'
#' @param g_init The prior distribution. Usually this is left unspecified and
#'   estimated from the data. However, it can be used in conjuction with
#'   \code{fix_g = TRUE} to fix the prior (useful, for example, to do
#'   computations with the "true" \code{g}). If \code{g_init} is specified but
#'   \code{fix_g = FALSE}, \code{g_init} specifies the initial value of \code{g}
#'   used during optimization.
#'
#' @param fix_g If \code{TRUE}, fix \code{g} at the specified value instead of
#'   estimating it.
#'
#' @param output A character vector indicating which values are to be returned.
#'   Options include:
#'     \describe{
#'       \item{\code{"result"}}{Summary results (posterior first and second
#'         moments).}
#'       \item{\code{"fitted_g"}}{The fitted prior.}
#'       \item{\code{"loglik"}}{The optimal log likelihood attained.}
#'       \item{\code{"lfsr"}}{A vector of local false sign rates.}
#'       \item{\code{"post_sampler"}}{A function that can be used to produce
#'         samples from the posterior. It takes a single parameter
#'         \code{nsamp}, the number of posterior samples to return per
#'         observation.}
#'      }
#'
#' @param prior_type The class of distributions from which the prior is to be
#'   estimated. See "Details" below.
#'
#' @param ... Additional parameters. All \code{ashr} prior types pass these
#'   parameters to \code{ashr::ash}.
#'
#' @examples
#' theta <- c(rep(0, 1000), rexp(1000)) # means
#' s <- rgamma(2000, 1, 1) # standard errors
#' x <- theta + rnorm(2000, 0, s) # observations
#' x.ebnm <- ebnm_point_normal(x, s)
#' pm <- x.ebnm$result$posterior_mean
#'
#' @export
#'
ebnm <- function(x,
                 s = 1,
                 mode = 0,
                 scale = "estimate",
                 g_init = NULL,
                 fix_g = FALSE,
                 output = output_default(),
                 prior_type = c("point_normal",
                                "point_laplace",
                                "normal",
                                "normal_scale_mixture",
                                "unimodal",
                                "unimodal_symmetric",
                                "unimodal_nonnegative",
                                "unimodal_nonpositive",
                                "ash"),
                 ...) {
  prior_type <- match.arg(prior_type)
  check_args(x, s, g_init, fix_g, output)

  if (prior_type == "point_normal") {
    retlist <- ebnm_point_normal(x, s, mode, scale, g_init, fix_g, output, ...)
  } else if (prior_type == "point_laplace") {
    retlist <- ebnm_point_laplace(x, s, mode, scale, g_init, fix_g, output, ...)
  } else if (prior_type == "normal") {
    retlist <- ebnm_normal(x, s, mode, scale, g_init, fix_g, output, ...)
  } else if (prior_type == "normal_scale_mixture") {
    retlist <- ebnm_ash_workhorse(x, s, mode, scale, g_init, fix_g, output,
                                  call = match.call(),
                                  mixcompdist = "normal", ...)
  } else if (prior_type == "unimodal") {
    retlist <- ebnm_ash_workhorse(x, s, mode, scale, g_init, fix_g, output,
                                  call = match.call(),
                                  mixcompdist = "halfuniform", ...)
  } else if (prior_type == "unimodal_symmetric") {
    retlist <- ebnm_ash_workhorse(x, s, mode, scale, g_init, fix_g, output,
                                  mixcompdist = "uniform",
                                  call = match.call(), ...)
  } else if (prior_type == "unimodal_nonnegative") {
    retlist <- ebnm_ash_workhorse(x, s, mode, scale, g_init, fix_g, output,
                                  call = match.call(),
                                  mixcompdist = "+uniform", ...)
  } else if (prior_type == "unimodal_nonpositive") {
    retlist <- ebnm_ash_workhorse(x, s, mode, scale, g_init, fix_g, output,
                                  call = match.call(),
                                  mixcompdist = "-uniform", ...)
  } else if (prior_type == "ash") {
    retlist <- ebnm_ash_workhorse(x, s, mode, scale, g_init, fix_g, output,
                                  call = match.call(), ...)
  }

  return(retlist)
}

check_args <- function(x, s, g_init, fix_g, output) {
  if (!(length(s) %in% c(1, length(x)))) {
    stop("Argument 's' must have either length 1 or the same length as ",
         "argument 'x'.")
  }

  if (all(is.infinite(s))) {
    stop("Standard errors cannot all be infinite.")
  }

  if (fix_g && is.null(g_init)) {
    stop("If g is fixed, then an initial g must be provided.")
  }

  if (!all(output %in% output_all())) {
    stop("Invalid argument to output. See function output_all() for a list ",
         "of valid outputs.")
  }
}
