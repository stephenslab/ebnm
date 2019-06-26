#' Solve the EBNM problem
#'
#' Solves the empirical Bayes normal means problem using a specified family of
#'   priors.
#'
#' @details TODO: update me.
#'
#' Given vectors of data \code{x} and standard errors \code{s},
#'   solve the "empirical Bayes normal means" (EBNM) problem, for various
#'   choices of prior family.
#'   The model is \deqn{x_j | \theta_j, s_j \sim N(\theta_j, s_j^2),}  and
#'   \deqn{\theta_j | s_j \sim g \in G}, where the distribution \eqn{g} is to be estimated.
#'   The distribution \eqn{g} is often referred to as the "prior distribution" for \eqn{\theta_j}
#'   and \eqn{G} is a specified family of prior distributions (several options for \eqn{G} are implemented, some
#'   parametric and others non-parametric;  see below for examples).
#'
#'   Solving the EBNM problem involves
#'   two steps. First, estimate \eqn{g \in  G}, by maximum marginal likelihood, yielding
#'   an estimate \deqn{\hat{g}:= \arg\max_{g \in G} L(g)} where
#'   \deqn{L(g):= \prod_j
#'   \int p(x_j | \theta_j, s_j)  g(d\theta_j);}
#'   Second, compute the posterior distributions \eqn{p(\theta_j | x_j, s_j, \hat{g})}, and/or summaries
#'   such as the posterior means and posterior second moments, etc.
#'
#' @param x A vector of observations.
#'
#' @param s A vector of standard deviations (or a scalar if all are equal).
#'
#' @param mode Scalar specifying the mode of the prior, g. Set to \code{"estimate"}
#'   to estimate it from the data.
#'
#' @param scale Scalar or vector, specifying the scale parameter(s) of the prior. The
#' precise interpretation of \code{scale} depends on \code{prior_type}. For \code{prior_type=normal, point_normal}
#' it is a scalar specifying the standard deviation of the normal component;
#' for \code{prior_type=point_laplace} it is a scalar specifying the rate parameter of the
#' Laplace component; for other prior types, which are implemented using the \code{\link[ashr]{ash}} function in the
#' \code{ashr} package,
#' it is a vector specifying the parameter \code{mixsd} to be passed to \code{\link[ashr]{ash}}.
#' Set to \code{"estimate"} to estimate it from the data (or to use the default \code{mixsd} in \code{\link[ashr]{ash}}).
#'
#' @param g_init The prior distribution, \eqn{g}. Usually this is left unspecified (NULL) and
#'   estimated from the data. However, it can be used in conjuction with
#'   \code{fix_g = TRUE} to fix the prior (useful, for example, to do
#'   computations with the "true" \code{g} in simulations). If \code{g_init} is specified but
#'   \code{fix_g = FALSE}, \code{g_init} specifies the initial value of \code{g}
#'   used during optimization. This has the side effect of fixing the
#'   \code{mode} and \code{scale} parameters for adaptive shrinkage (ash) priors.
#'
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g}=\code{g_init} instead of
#'   estimating it.
#'
#' @param output A character vector indicating which values are to be returned.
#'   Options include:
#'     \describe{
#'       \item{\code{"result"}}{Summary results (posterior first and second
#'         moments).}
#'       \item{\code{"fitted_g"}}{The fitted prior \eqn{\hat{g}}.}
#'       \item{\code{"loglik"}}{The optimal log likelihood attained, \eqn{L(\hat{g}}.}
#'       \item{\code{"lfsr"}}{A vector of local false sign rates.}
#'       \item{\code{"post_sampler"}}{A function that can be used to produce
#'         samples from the posterior. It takes a single parameter
#'         \code{nsamp}, the number of posterior samples to return per
#'         observation.}
#'      }
#'
#' @param control A list of control parameters to be passed to the optimization
#'   function (\code{nlm} for normal, point-normal, and point-Laplace priors
#'   and, unless specified otherwise, \code{mixsqp::mixsqp} for ash priors).
#'
#' @param prior_type A character string that specifies the prior family \eqn{G}. See "Details" below.
#'
#' @param ... Additional parameters. \code{unimodal_} prior types pass these
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
                 control = NULL,
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

  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        control = control,
                        prior_type = prior_type,
                        call = match.call(),
                        ...))
}

# All exported functions call into ebnm_workhorse. This allows argument
#   checking to be done in a single place.
#
ebnm_workhorse <- function(x,
                           s,
                           mode,
                           scale,
                           g_init,
                           fix_g,
                           output,
                           control,
                           prior_type,
                           call,
                           ...) {
  check_args(x, s, g_init, fix_g, output)
  mode <- handle_mode_parameter(mode)
  scale <- handle_scale_parameter(scale)
  if (is.null(control)) {
    control <- list()
  }

  if (prior_type == "point_normal") {
    retlist <- ebnm_pn_workhorse(x = x,
                                 s = s,
                                 mode = mode,
                                 scale = scale,
                                 g_init = g_init,
                                 fix_g = fix_g,
                                 output = output,
                                 control = control,
                                 pointmass = TRUE,
                                 call = call)
  } else if (prior_type == "point_laplace") {
    retlist <- ebnm_pl_workhorse(x = x,
                                 s = s,
                                 mode = mode,
                                 scale = scale,
                                 g_init = g_init,
                                 fix_g = fix_g,
                                 output = output,
                                 control = control,
                                 call = call)
  } else if (prior_type == "normal") {
    retlist <- ebnm_pn_workhorse(x = x,
                                 s = s,
                                 mode = mode,
                                 scale = scale,
                                 g_init = g_init,
                                 fix_g = fix_g,
                                 output = output,
                                 control = control,
                                 pointmass = FALSE,
                                 call = call)
  } else if (prior_type == "normal_scale_mixture") {
    retlist <- ebnm_normal_mix_workhorse(x = x,
                                         s = s,
                                         mode = mode,
                                         scale = scale,
                                         g_init = g_init,
                                         fix_g = fix_g,
                                         output = output,
                                         control = control,
                                         pointmass = TRUE,
                                         grid_mult = sqrt(2),
                                         call = call)
  } else if (prior_type == "unimodal") {
    retlist <- ebnm_ash_workhorse(x = x,
                                  s = s,
                                  mode = mode,
                                  scale = scale,
                                  g_init = g_init,
                                  fix_g = fix_g,
                                  output = output,
                                  control = control,
                                  call = call,
                                  mixcompdist = "halfuniform",
                                  ...)
  } else if (prior_type == "unimodal_symmetric") {
    retlist <- ebnm_ash_workhorse(x = x,
                                  s = s,
                                  mode = mode,
                                  scale = scale,
                                  g_init = g_init,
                                  fix_g = fix_g,
                                  output = output,
                                  control = control,
                                  call = call,
                                  mixcompdist = "uniform",
                                  ...)
  } else if (prior_type == "unimodal_nonnegative") {
    retlist <- ebnm_ash_workhorse(x = x,
                                  s = s,
                                  mode = mode,
                                  scale = scale,
                                  g_init = g_init,
                                  fix_g = fix_g,
                                  output = output,
                                  control = control,
                                  call = call,
                                  mixcompdist = "+uniform",
                                  ...)
  } else if (prior_type == "unimodal_nonpositive") {
    retlist <- ebnm_ash_workhorse(x = x,
                                  s = s,
                                  mode = mode,
                                  scale = scale,
                                  g_init = g_init,
                                  fix_g = fix_g,
                                  output = output,
                                  control = control,
                                  call = call,
                                  mixcompdist = "-uniform",
                                  ...)
  } else if (prior_type == "ash") {
    retlist <- ebnm_ash_workhorse(x = x,
                                  s = s,
                                  mode = mode,
                                  scale = scale,
                                  g_init = g_init,
                                  fix_g = fix_g,
                                  output = output,
                                  control = control,
                                  call = call,
                                  ...)
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

handle_mode_parameter <- function(mode) {
  # Allow partial matching.
  if (identical(pmatch(mode, "estimate"), 1L)) {
    mode <- "estimate"
  } else if (!(is.numeric(mode) && length(mode) == 1 && is.finite(mode))) {
    stop("Argument 'mode' must be either 'estimate' or a numeric value.")
  }
  return(mode)
}

handle_scale_parameter <- function(scale) {
  # Allow partial matching.
  if (identical(pmatch(scale, "estimate"), 1L)) {
    scale <- "estimate"
  }
  else if (!(is.numeric(scale) && all(is.finite(scale)))) {
    stop("Argument 'scale' must be either 'estimate' or numeric.")
  }
  return(scale)
}
