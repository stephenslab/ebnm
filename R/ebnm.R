#' Solve the EBNM problem
#'
#' Solves the Empirical Bayes Normal Means problem using either a point-normal
#'   or a point-laplace prior.
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
#'   }
#'
#' @param g The prior distribution (a list with elements \code{pi0}, \code{a},
#'   and \code{mu}). Usually this is left unspecified and estimated from the
#'   data. However, it can be used in conjuction with \code{fixg = TRUE} to
#'   fix the prior (useful, for example, to do computations with the "true"
#'   \code{g}). If \code{g} is specified but \code{fixg = FALSE}, \code{g}
#'   specifies the initial value of \code{g} used during optimization.
#'
#' @param fixg If \code{TRUE}, fix \code{g} at the specified value instead of
#'   estimating it. This overrides any settings of parameters \code{fix_pi0}
#'   and \code{fix_mu}.
#'
#' @param fix_pi0 If \code{TRUE}, fix \code{g$pi0} at the specified value.
#'   \code{g$a} will be estimated from the data and \code{g$mu} will be fixed
#'   or estimated according to the setting of parameter \code{fix_mu}. This
#'   option has not yet been implemented for the point-laplace prior.
#'
#' @param fix_mu If \code{TRUE}, fix \code{g$mu} at the specified value (or
#'   at zero if \code{g$mu} has not been specified). \code{g$a} will be
#'   estimated and \code{g$mu} will be fixed or estimated depending on the
#'   setting of parameter \code{fix_pi0}. This option is only available for
#'   the point-normal prior.
#'
#' @param norm The normalization factor to divide \code{x} and \code{s}
#'   by before running optimization. This should not affect results, but
#'   it can improve numerical stability when \code{x} and \code{s} are very
#'   small.
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
#'      }
#'
#' @examples
#' theta = c(rep(0, 1000), rexp(1000)) # means
#' s = rgamma(2000, 1, 1) # standard errors
#' x = theta + rnorm(2000, 0, s) # observations
#' x.ebnm = ebnm(x, s, "point_normal")
#' ashr::get_pm(x.ebnm) # posterior mean
#'
#' @export
#'
ebnm <- function(x,
                 s = 1,
                 prior_type = c("point_normal", "point_laplace"),
                 g = list(),
                 fixg = FALSE,
                 fix_pi0 = FALSE,
                 fix_mu = TRUE,
                 norm = NULL,
                 control = NULL,
                 output = NULL) {
  prior_type <- match.arg(prior_type)

  if (prior_type == "point_normal") {
    retlist <- ebnm_point_normal(x, s, g, fixg,
                                 fix_pi0, fix_mu, norm, control, output)
  } else {
    if (!fix_mu || (!is.null(g$mu) && g$mu != 0)) {
      stop("Currently, 'mu' must be fixed at zero for 'point_laplace' ",
           " priors.")
    }
    retlist <- ebnm_point_laplace(x, s, g, fixg, output)
  }

  return(retlist)
}

# Argument checks shared by main exported functions (ebnm, ebnm_point_normal,
#   and ebnm_point_laplace).
check_args <- function(x, s, g, fixg, output) {
  if (!(length(s) %in% c(1, length(x)))) {
    stop("Argument 's' must have either length 1 or the same length as ",
         "argument 'x'.")
  }
  if (all(is.infinite(s))) {
    stop("Standard errors cannot all be infinite.")
  }

  if (fixg && (is.null(g$a) || is.null(g$pi0))) {
    stop("Must specify g$pi0 and g$a if fixg = TRUE.")
  }
  if (!is.null(g$a) && (g$a <= 0)) {
    stop("Invalid choice of g$a.")
  }
  if (!is.null(g$pi0) && (g$pi0 < 0 || g$pi0 > 1)) {
    stop("Invalid choice of g$pi0.")
  }

  if (!all(output %in% c("result", "fitted_g", "loglik", "post_sampler"))) {
    stop("Argument 'output' must consist of 'result', 'fitted_g', 'loglik', ",
         "and/or 'post_sampler'.")
  }
}
