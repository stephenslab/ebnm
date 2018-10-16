#' @title Solve the EBNM problem with normal prior
#'
#' @description This function solves the Empirical Bayes Normal Means
#'   problem with a normal prior.
#'
#' @details Given vectors of data \code{x} and standard errors \code{s},
#'   solve the EBNM problem with a normal prior. The model is
#'  \deqn{x_j \sim N(\theta_j, s_j^2),} where \eqn{s_j} are given and
#'   \eqn{\theta_j \sim N(\mu, 1/a)}. \eqn{\mu} and \eqn{a} are estimated by
#'   marginal maximum likelihood.
#'
#' @param x A vector of observations.
#'
#' @param s A vector of standard deviations (or a scalar if all are
#'   equal).
#'
#' @param g The prior distribution (a list with elements \code{mu} and
#'   \code{a}). Usually this is unspecified (\code{NULL}) and estimated
#'   from the data. However, it can be used in conjuction with
#'   \code{fixg = TRUE} to specify the prior to use (useful, for example,
#'   to do computations with the "true" \code{g}). Or, if \code{g} is
#'   specified but \code{fixg = FALSE}, \code{g} specifies the initial
#'   value of \code{g} used before optimization.
#'
#' @param fixg If \code{TRUE}, use the specified \code{g} instead of
#'   estimating it.
#'
#' @param fix_mu If \code{TRUE}, \code{g$mu} is fixed at the supplied
#'   value and \code{g$a} is estimated from the data. \code{fixg = TRUE}
#'   overrides \code{fix_mu = TRUE}. That is, if both are \code{TRUE}
#'   then both \code{mu} and \code{a} are fixed at the supplied values.
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
#'   \item{"fitted_g"}{The fitted prior (a list with elements \code{mu}
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
#' theta = rnorm(2000, 3, 1) # means
#' s = rgamma(2000, 1, 1) # standard errors
#' x = theta + rnorm(2000, 0, s) # observations
#' x.ebnm = ebnm_normal(x, s)
#' ashr::get_pm(x.ebnm) # posterior mean
#'
#' @export
#'
ebnm_normal <- function (x,
                         s = 1,
                         g = NULL,
                         fixg = FALSE,
                         fix_mu = FALSE,
                         norm = NULL,
                         control = NULL,
                         output = c("result", "fitted_g", "loglik")) {
  output <- set_output(output)
  
  if (fixg && is.null(g)) {
    stop("Must specify g if fixg = TRUE")
  }
  if (fix_mu && is.null(g$mu)) {
    stop("Must specify g$mu if fix_mu = TRUE")
  }
  if (!is.null(g$a) && (g$a <= 0)) {
    stop("Invalid choice of g$a")
  }
  if (!is.null(g$mu) && is.infinite(g$mu)) {
    stop("Invalid choice of g$mu")
  }
  
  if (all(is.infinite(s))) {
    stop("Impossible to fit g when all standard errors are infinite")
  }
  
  # Scale for stability, but need to be careful with log-likelihood:
  if (is.null(norm)) {
    pos_idx <- (is.finite(s) & s > 0)
    if (sum(pos_idx) == 0) {
      norm <- 1
    } else {
      norm <- mean(s[pos_idx])
    }
  }
  
  s <- s / norm
  x <- x / norm
  if (!is.null(g) && !is.null(g$a)) {
    g$a <- g$a * norm^2
  }
  if (!is.null(g) && !is.null(g$mu)) {
    g$mu <- g$mu / norm
  }
  
  # Don't use data with infinite SEs when estimating g:
  x_subset <- x[is.finite(s)]
  s_subset <- s[is.finite(s)]
  
  if (!fixg) {
    if (fix_mu) {
      g <- mle_normal_logscale_fixed_mu(x_subset, s_subset, g, control)
    } else {
      g <- mle_normal_logscale_grad(x_subset, s_subset, g, control)
    }
  }
  
  mu <- g$mu
  a <- g$a
  
  # Compute return values, taking care to adjust back to original scale:
  retlist <- list()
  
  if ("result" %in% output) {
    result <- compute_summary_results_normal(x, s, mu, a)
    result$PosteriorMean <- result$PosteriorMean * norm
    result$PosteriorMean2 <- result$PosteriorMean2 * norm^2
    retlist <- c(retlist, list(result = result))
  }
  
  if ("fitted_g" %in% output) {
    fitted_g <- list(mu = g$mu * norm, a = g$a / norm^2)
    retlist <- c(retlist, list(fitted_g = fitted_g))
  }
  
  if ("loglik" %in% output) {
    loglik <- loglik_normal(x_subset, s_subset, mu, a)
    loglik <- loglik - length(x_subset) * log(norm)
    retlist <- c(retlist, list(loglik = loglik))
  }
  
  if ("post_sampler" %in% output) {
    retlist <- c(retlist, list(post_sampler = function(nsamp) {
      post_sampler_normal(x, s, mu, a, nsamp) * norm
    }))
  }
  
  return(retlist)
}
