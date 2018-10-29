#' @title Solve the EBNM problem with point-normal or point-laplace prior
#'
#' @description This function solves the Empirical Bayes Normal Means
#'   problem with a point-normal or point-laplace prior.
#'
#' @details Given vectors of data \code{x} and standard errors \code{s},
#'   solve the EBNM problem with a "point-normal" or "point-laplace" prior. The model is
#'  \deqn{x_j \sim N(\theta_j, s_j^2),} where \eqn{s_j} are given and
#'   \eqn{\theta_j \sim g}, with \eqn{g}, either a mixture of a point mass at
#'   \eqn{\mu} and a normal distribution: \deqn{\theta_j \sim \pi_0 \delta_\mu
#'   + (1 - \pi_0)N(\mu, 1/a)}, or a mixture of a point mass at 0 and a laplace distribution:
#'   \deqn{\theta_j \sim \pi_0 \delta_0 + (1 - \pi_0)DExp(a)}. \eqn{\pi_0} and \eqn{a} 
#'   (and \eqn{\mu}, for the nornal case) are estimated by marginal maximum likelihood.
#'
#' @param x A vector of observations.
#'
#' @param s A vector of standard deviations (or a scalar if all are
#'   equal).
#'   
#' @param prior_type The specified prior for the non-null component of \eqn{\theta}.
#' Options include:
#' \describe{
#'   \item{"point_normal"}{Prior is a mixture of a point-mass and a normal distribution.}
#'   \item{"point_laplace"}{Prior is a mixture of a point-mass and a laplace distribution.}
#' }
#'
#' @param g The prior distribution (a list with elements \code{pi0}, 
#'   \code{a}, and \code{mu}). Usually this is unspecified (\code{NULL}) and estimated
#'   from the data. However, it can be used in conjuction with
#'   \code{fixg = TRUE} to specify the prior to use (useful, for example,
#'   to do computations with the "true" \code{g}). Or, if \code{g} is
#'   specified but \code{fixg = FALSE}, \code{g} specifies the initial
#'   value of \code{g} used before optimization.
#'
#' @param fixg If \code{TRUE}, use the specified \code{g} instead of
#'   estimating it.
#'
#' @param fix_pi0 If \code{TRUE}, \code{g$pi0} is fixed at the supplied
#'   value and \code{g$a} and \code{g$mu} are estimated from the data. 
#'   \code{fixg = TRUE} overrides \code{fix_pi0 = TRUE}. That is, if both 
#'   are \code{TRUE} then all of \code{pi0} and \code{a} and \code{mu} are fixed 
#'   at the supplied values.
#'   
#' @param fix_mu If \code{TRUE}, \code{g$mu} is fixed at the supplied
#'   value and \code{g$a} and \code{g$pi0} is estimated from the data. 
#'   \code{fixg = TRUE} overrides \code{fix_mu = TRUE}. That is, if both are 
#'   \code{TRUE} then all of \code{pi0} and \code{a} and \code{mu} fixed at 
#'   the supplied values.
#'   At this time, this option is only available for the point-normal prior.
#'   If one wishes to fix the mean at \code{mu}, one can simply subtract 
#'   \code{mu} from all observations \code{x}.
#'
#' @param norm The normalization factor to divide \code{x} and \code{s}
#'   by before running optimization (this should not affect results, but
#'   it improves numerical stability when \code{x} and \code{s} are
#'   tiny).
#'
#' @param control A list of control parameters to be passed to
#'   \code{optim}.
#'
#' @param output A vector of strings indicating which values are to be
#' returned. Options include:
#' \describe{
#'   \item{"result"}{Summary results (posterior means \eqn{E \theta_j}
#'     and posterior values of \eqn{E \theta_j^2}).}
#'   \item{"fitted_g"}{The fitted prior (a list with elements \code{pi0}
#'     and \code{a}).}
#'   \item{"loglik"}{The optimal log likelihood attained.}
#'   \item{"post_sampler"}{A function that can be used to produce
#'     samples from the posterior. It takes a single parameter
#'     \code{nsamp}, the number of posterior samples to return per
#'     observation.}
#'  }
#'
#' @return A list with elements specified by parameter \code{output}.
#'
#' @examples
#' mu = c(rep(0, 1000), rexp(1000)) # means
#' s = rgamma(2000, 1, 1) # standard errors
#' x = mu + rnorm(2000, 0, s) # observations
#' x.ebnm = ebnm(x, s, "point_normal")
#' ashr::get_pm(x.ebnm) # posterior mean
#'
#' @export
ebnm <- function(x, 
                 s = 1, 
                 prior_type,
                 g = list(mu = 0),
                 fixg = FALSE,
                 fix_pi0 = FALSE,
                 fix_mu = TRUE,
                 norm = NULL,
                 control = NULL,
                 output = c("result", "fitted_g", "loglik")) {
  
  output <- set_output(output)
  
  # make sure s is either of length 1 or the same length of x
  if (!(length(s) %in% c(1, length(x)))) {
    stop("Argument 's' must have either length 1, or the same length of argument 'x'.")
  }
  
  # make sure prior_type and output are correctly specified
  if (!(prior_type %in% c("point_normal", "point_laplace"))) {
    stop("Argument 'prior_type' must be one of 'point_normal' or 'point_laplace")
  }
  if (any(!(output %in% c("result", "fitted_g", "loglik", "post_sampler")))) {
    stop("Argument 'output' must consist of 'result', 'fitted_g', 'loglik', and/or 'post_sampler'")
  }
  
  # currently, estimating mu not supported for point_laplace
  if ((prior_type == "point_laplace") && !fix_mu) {
    warning("Prior 'point_laplace' currently does not support estimating 'mu'. Fixing 'mu' to be 0")
  }
  
  if (prior_type == "point_normal") {
    retlist = ebnm_point_normal(x, s, g, fixg, fix_pi0, fix_mu, norm, control, output)
  } else {
    if (!is.null(g)) {
      if (names(g) == "mu") { # if only mu supplied
        g = NULL
      } else{
        g$mu = NULL
      }
    }
    retlist = ebnm_point_laplace(x, s, g, fixg, output)
    if ("fitted_g" %in% output) {
      retlist$fitted_g$mu = 0
    }
  }
  
  return(retlist)
  
}

