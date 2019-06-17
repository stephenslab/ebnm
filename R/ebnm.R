#' Solve the EBNM problem
#'
#' Solves the Empirical Bayes Normal Means problem using a point-normal,
#'   point-laplace, or normal prior.
#'
#' @details Given vectors of data \code{x} and standard errors \code{s},
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
#' @param prior_type The type of prior to estimate. Options include:
#'   \describe{
#'     \item{\code{"point_normal"}}{Prior is a mixture of a point mass and a
#'       normal distribution.}
#'     \item{\code{"point_laplace"}}{Prior is a mixture of a point mass and a
#'       laplace distribution.}
#'     \item{\code{"normal"}}{Prior is a normal distribution (with no point
#'       mass).}
#'   }
#'
#' @param g The prior distribution (a list with elements \code{pi0}, \code{a},
#'   and \code{mu}). Usually this is left unspecified and estimated from the
#'   data. However, it can be used in conjuction with \code{fix_g = TRUE} to
#'   fix the prior (useful, for example, to do computations with the "true"
#'   \code{g}). If \code{g} is specified but \code{fix_g = FALSE}, \code{g}
#'   specifies the initial value of \code{g} used during optimization.
#'
#' @param fix_g If \code{TRUE}, fix \code{g} at the specified value instead of
#'   estimating it. This overrides any settings of parameters \code{fix_pi0},
#'   \code{fix_a}, and \code{fix_mu}.
#'
#' @param fix_pi0,fix_a,fix_mu For a point-normal prior, any combination of
#'   \code{pi0}, \code{a}, and \code{mu} can be fixed. This functionality has
#'   not yet been implemented for the point-laplace prior.
#'
#' @param control A list of control parameters to be passed to \code{optim}.
#'
#' @param output A character vector indicating which values are to be returned.
#'   Default values are set by function \code{ebnm:::set_output}. Options
#'   include:
#'     \describe{
#'       \item{\code{"result"}}{Summary results (posterior means
#'         \eqn{E \theta_j} and posterior values of \eqn{E \theta_j^2}).}
#'       \item{\code{"fitted_g"}}{The fitted prior (a list with elements
#'         \code{pi0}, \code{a}, and \code{mu}).}
#'       \item{\code{"loglik"}}{The optimal log likelihood attained.}
#'       \item{\code{"post_sampler"}}{A function that can be used to produce
#'         samples from the posterior. It takes a single parameter
#'         \code{nsamp}, the number of posterior samples to return per
#'         observation.}
#'       \item{\code{"lfsr"}}{A vector of local false sign rates.}
#'      }
#'
#' @examples
#' theta <- c(rep(0, 1000), rexp(1000)) # means
#' s <- rgamma(2000, 1, 1) # standard errors
#' x <- theta + rnorm(2000, 0, s) # observations
#' x.ebnm <- ebnm(x, s, "point_normal")
#' pm <- x.ebnm$PosteriorMean
#'
#' @export
#'
ebnm <- function(x,
                 s = 1,
                 g_init = NULL,
                 fix_g = FALSE,
                 output = output_default(),
                 prior_type = c("point_normal",
                                "normal",
                                "point_laplace",
                                "ash"),
                 ...) {
  prior_type <- match.arg(prior_type)
  check_args(x, s, g_init, fix_g, output)

  if (prior_type == "point_normal") {
    retlist <- ebnm_point_normal(x, s, g_init, fix_g, output, ...)
  } else if (prior_type == "normal") {
    retlist <- ebnm_normal(x, s, g_init, fix_g, output, ...)
  } else if (prior_type == "point_laplace") {
    retlist <- ebnm_point_laplace(x, s, g_init, fix_g, output, ...)
  } else if (prior_type == "ash") {
    retlist <- ebnm_ash(x, s, g_init, fix_g, output, ...)
  }

  return(retlist)
}

output_default <- function() {
  return(c("result", "fitted_g", "loglik"))
}

#' @export
output_all <- function() {
  return(c("result", "fitted_g", "loglik", "post_sampler", "lfsr"))
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
