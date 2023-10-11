#' Solve the EBNM problem using point-normal priors
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   point-normal priors (the family of mixtures where one component is a point
#'   mass at \eqn{\mu} and the other is a normal distribution centered at
#'   \eqn{\mu}). Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "point_normal"}. For details about the model, see
#'   \code{\link{ebnm}}.
#'
#' @param x A vector of observations. Missing observations (\code{NA}s) are
#'   not allowed.
#'
#' @param s A vector of standard errors (or a scalar if all are equal).
#'   Standard errors may not be exactly zero, and missing standard errors are
#'   not allowed.
#'
#' @param mode A scalar specifying the mode of the prior \eqn{g} or
#'   \code{"estimate"} if the mode is to be estimated from the data.
#'
#' @param scale A scalar specifying the standard deviation of the normal
#'   component or \code{"estimate"} if the standard deviation is to be estimated
#'   from the data.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link[ashr]{normalmix}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{normalmix}.
#'
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
#'   of estimating it.
#'
#' @param output A character vector indicating which values are to be returned.
#'   Function \code{ebnm_output_default()} provides the default return values, while
#'   \code{ebnm_output_all()} lists all possible return values. See \strong{Value}
#'   below.
#'
#' @param optmethod A string specifying which optimization function is to be
#'   used. Options include \code{"nlm"}, \code{"lbfgsb"} (which calls
#'   \code{optim} with \code{method = "L-BFGS-B"}), and \code{"trust"} (which
#'   calls into package \code{trust}). Other options are \code{"nohess_nlm"},
#'   \code{"nograd_nlm"}, and \code{"nograd_lbfgsb"}, which use numerical
#'   approximations rather than exact expressions for the Hessian and (for
#'   the latter two) the gradient. The default option is \code{"nohess_nlm"}.
#'
#' @param control A list of control parameters to be passed to the
#'   optimization function specified by parameter \code{optmethod}.
#'
#' @return An \code{ebnm} object. Depending on the argument to \code{output}, the
#'   object is a list containing elements:
#'     \describe{
#'       \item{\code{data}}{A data frame containing the observations \code{x}
#'         and standard errors \code{s}.}
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, standard deviations, second moments, and local false sign
#'         rates).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}}.}
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained,
#'         \eqn{L(\hat{g})}.}
#'       \item{\code{posterior_sampler}}{A function that can be used to
#'         produce samples from the posterior. The sampler takes a single
#'         parameter \code{nsamp}, the number of posterior samples to return per
#'         observation.}
#'      }
#'    S3 methods \code{coef}, \code{confint}, \code{fitted}, \code{logLik},
#'    \code{nobs}, \code{plot}, \code{predict}, \code{print}, \code{quantile},
#'    \code{residuals}, \code{simulate}, \code{summary}, and \code{vcov}
#'    have been implemented for \code{ebnm} objects. For details, see the
#'    respective help pages, linked below under \strong{See Also}.
#'
#' @seealso See \code{\link{ebnm}} for examples of usage and model details.
#'
#'   Available S3 methods include \code{\link{coef.ebnm}},
#'   \code{\link{confint.ebnm}},
#'   \code{\link{fitted.ebnm}}, \code{\link{logLik.ebnm}},
#'   \code{\link{nobs.ebnm}}, \code{\link{plot.ebnm}},
#'   \code{\link{predict.ebnm}}, \code{\link{print.ebnm}},
#'   \code{\link{print.summary.ebnm}}, \code{\link{quantile.ebnm}},
#'   \code{\link{residuals.ebnm}}, \code{\link{simulate.ebnm}},
#'   \code{\link{summary.ebnm}}, and \code{\link{vcov.ebnm}}.
#'
#' @export
#'
ebnm_point_normal <- function(x,
                              s = 1,
                              mode = 0,
                              scale = "estimate",
                              g_init = NULL,
                              fix_g = FALSE,
                              output = ebnm_output_default(),
                              optmethod = NULL,
                              control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = optmethod,
                        control = control,
                        prior_family = "point_normal",
                        call = match.call()))
}

#' Solve the EBNM problem using point-Laplace priors
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   point-Laplace priors (the family of mixtures where one component is a point
#'   mass at \eqn{\mu} and the other is a double-exponential distribution
#'   centered at \eqn{\mu}). Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "point_laplace"}. For details about the model, see
#'   \code{\link{ebnm}}.
#'
#' @inherit ebnm_point_normal
#'
#' @param scale A scalar specifying the scale parameter of the Laplace
#'   component or \code{"estimate"} if the scale is to be estimated
#'   from the data.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link{laplacemix}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{laplacemix}.
#'
#' @export
#'
ebnm_point_laplace <- function(x,
                               s = 1,
                               mode = 0,
                               scale = "estimate",
                               g_init = NULL,
                               fix_g = FALSE,
                               output = ebnm_output_default(),
                               optmethod = NULL,
                               control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = optmethod,
                        control = control,
                        prior_family = "point_laplace",
                        call = match.call()))
}

#' Solve the EBNM problem using point-exponential priors
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   point-exponential priors (the family of mixtures where one component is a
#'   point mass at \eqn{\mu} and the other is a (nonnegative) exponential
#'   distribution with mode \eqn{\mu}). Identical to function \code{\link{ebnm}}
#'   with argument \code{prior_family = "point_exponential"}. For details about
#'   the model, see \code{\link{ebnm}}.
#'
#' @inherit ebnm_point_normal
#'
#' @param scale A scalar specifying the scale parameter of the exponential
#'   component or \code{"estimate"} if the scale is to be estimated
#'   from the data.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link{gammamix}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{gammamix}.
#'
#' @export
#'
ebnm_point_exponential <- function(x,
                                   s = 1,
                                   mode = 0,
                                   scale = "estimate",
                                   g_init = NULL,
                                   fix_g = FALSE,
                                   output = ebnm_output_default(),
                                   optmethod = NULL,
                                   control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = optmethod,
                        control = control,
                        prior_family = "point_exponential",
                        call = match.call()))
}

#' Solve the EBNM problem using normal priors
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   normal distributions. Identical to function \code{\link{ebnm}} with
#'   argument \code{prior_family = "normal"}. For details about the model, see
#'   \code{\link{ebnm}}.
#'
#' @inherit ebnm_point_normal
#'
#' @param scale A scalar specifying the standard deviation of the normal prior
#'   or \code{"estimate"} if the standard deviation is to be estimated from
#'   the data.
#'
#' @export
#'
ebnm_normal <- function(x,
                        s = 1,
                        mode = 0,
                        scale = "estimate",
                        g_init = NULL,
                        fix_g = FALSE,
                        output = ebnm_output_default(),
                        optmethod = NULL,
                        control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = optmethod,
                        control = control,
                        prior_family = "normal",
                        call = match.call()))
}

#' Solve the EBNM problem using horseshoe priors
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   \link{horseshoe} distributions. Identical to function \code{\link{ebnm}}
#'   with argument \code{prior_family = "horseshoe"}. For details about the
#'   model, see \code{\link{ebnm}}.
#'
#' @inherit ebnm_point_normal
#'
#' @param s A \emph{scalar} specifying the standard error of the observations
#'   (observations must be homoskedastic).
#'
#' @param scale A scalar corresponding to \eqn{s\tau} in the usual
#'   parametrization of the \code{\link{horseshoe}} distribution, or
#'   \code{"estimate"} if this parameter is to be estimated from the data.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link{horseshoe}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{horseshoe}.
#'
#' @param control A list of control parameters to be passed to function
#'   \code{\link[stats]{optimize}}.
#'
#' @return An \code{ebnm} object. Depending on the argument to \code{output}, the
#'   object is a list containing elements:
#'     \describe{
#'       \item{\code{data}}{A data frame containing the observations \code{x}
#'         and standard errors \code{s}.}
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, standard deviations, second moments, and local false sign
#'         rates).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}}.}
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained,
#'         \eqn{L(\hat{g})}.}
#'       \item{\code{posterior_sampler}}{A function that can be used to
#'         produce samples from the posterior. The function takes parameters
#'         \code{nsamp}, the number of posterior samples to return per
#'         observation, and \code{burn}, the number of burn-in samples to
#'         discard (an MCMC sampler is used).}
#'      }
#'    S3 methods \code{coef}, \code{confint}, \code{fitted}, \code{logLik},
#'    \code{nobs}, \code{plot}, \code{predict}, \code{print}, \code{quantile},
#'    \code{residuals}, \code{simulate}, \code{summary}, and \code{vcov}
#'    have been implemented for \code{ebnm} objects. For details, see the
#'    respective help pages, linked below under \strong{See Also}.
#'
#' @export
#'
ebnm_horseshoe <- function(x,
                           s = 1,
                           scale = "estimate",
                           g_init = NULL,
                           fix_g = FALSE,
                           output = ebnm_output_default(),
                           control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = 0,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "horseshoe",
                        call = match.call()))
}

#' Solve the EBNM problem using scale mixtures of normals
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   scale mixtures of normals. Identical to function \code{\link{ebnm}}
#'   with argument \code{prior_family = "normal_scale_mixture"}. For details
#'   about the model, see \code{\link{ebnm}}.
#'
#' @inherit ebnm_point_normal
#'
#' @param scale The nonparametric family of scale mixtures of normals is
#'   approximated via a finite mixture of normal distributions
#'   \deqn{\pi_1 N(\mu, \sigma_1^2) + \ldots + \pi_K N(\mu, \sigma_K^2),}
#'   where parameters \eqn{\pi_k} are estimated and the grid of standard
#'   deviations \eqn{(\sigma_1, \ldots, \sigma_K)} is fixed in advance. By
#'   making the grid sufficiently dense, one can obtain an arbitrarily good
#'   approximation. The grid can be specified by the user via parameter
#'   \code{scale}, in which case the argument should be the vector of
#'   standard deviations \eqn{(\sigma_1, \ldots, \sigma_K)}; alternatively,
#'   if \code{scale = "estimate"}, then
#'   \code{ebnm} sets the grid via function \code{\link{ebnm_scale_normalmix}}.
#'   Note that \code{ebnm} sets the grid differently from
#'   function \code{\link[ashr]{ash}}. To use the \code{ash} grid, set
#'   \code{scale = "estimate"} and pass in \code{gridmult} as an additional
#'   parameter. See \code{\link[ashr]{ash}} for defaults and details.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. This has
#'   the side effect of fixing the \code{mode} and \code{scale} parameters. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link[ashr]{normalmix}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{normalmix}.
#'
#' @param control A list of control parameters to be passed to optimization
#'   function \code{\link[mixsqp]{mixsqp}}.
#'
#' @param ... When parameter \code{gridmult} is set, an
#'   \code{\link[ashr]{ash}}-style grid will be used instead of the default
#'   \code{ebnm} grid (see parameter \code{scale} above). Other additional
#'   parameters are ignored.
#'
#' @export
#'
ebnm_normal_scale_mixture <- function(x,
                                      s = 1,
                                      mode = 0,
                                      scale = "estimate",
                                      g_init = NULL,
                                      fix_g = FALSE,
                                      output = ebnm_output_default(),
                                      control = NULL,
                                      ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "normal_scale_mixture",
                        call = match.call(),
                        ...))
}

#' Solve the EBNM problem using unimodal distributions
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   all unimodal distributions. Identical to function \code{\link{ebnm}}
#'   with argument \code{prior_family = "unimodal"}. For details
#'   about the model, see \code{\link{ebnm}}.
#'
#' @inherit ebnm_normal_scale_mixture
#'
#' @param scale The nonparametric family of unimodal distributions is
#'   approximated via a finite mixture of uniform distributions
#'   \deqn{\pi_1^l \mathrm{Unif}(\mu - a_1, \mu) + \pi_1^u \mathrm{Unif}(\mu, \mu + a_1)
#'   + \ldots + \pi_K^l \mathrm{Unif}(\mu - a_K, \mu) + \pi_K^u \mathrm{Unif}(\mu, \mu + a_K),}
#'   where parameters \eqn{\pi_k^l} and \eqn{\pi_k^u} are estimated and the grid
#'   of lengths \eqn{(a_1, \ldots, a_K)} is fixed in advance. By
#'   making the grid sufficiently dense, one can obtain an arbitrarily good
#'   approximation. The grid can be specified by the user via parameter
#'   \code{scale}, in which case the argument should be the vector of
#'   lengths \eqn{(a_1, \ldots, a_K)}; alternatively, if
#'   \code{scale = "estimate"}, then \code{ebnm} sets the grid via function
#'   \code{\link{ebnm_scale_unimix}}.
#'   Note that \code{ebnm} sets the grid differently from
#'   function \code{\link[ashr]{ash}}. To use the \code{ash} grid, set
#'   \code{scale = "estimate"} and pass in \code{gridmult} as an additional
#'   parameter. See \code{\link[ashr]{ash}} for defaults and details.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. This has
#'   the side effect of fixing the \code{mode} and \code{scale} parameters. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link[ashr]{unimix}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{unimix}.
#'
#' @param ... Additional parameters to be passed to function
#'   \code{\link[ashr]{ash}} in package \code{ashr}.
#'
#' @export
#'
ebnm_unimodal <- function(x,
                          s = 1,
                          mode = 0,
                          scale = "estimate",
                          g_init = NULL,
                          fix_g = FALSE,
                          output = ebnm_output_default(),
                          control = NULL,
                          ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "unimodal",
                        call = match.call(),
                        ...))
}

#' Solve the EBNM problem using symmetric unimodal distributions
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   symmetric unimodal distributions. Identical to function \code{\link{ebnm}}
#'   with argument \code{prior_family = "unimodal_symmetric"}. For details
#'   about the model, see \code{\link{ebnm}}.
#'
#' @inherit ebnm_unimodal
#'
#' @param scale The nonparametric family of symmetric unimodal distributions is
#'   approximated via a finite mixture of uniform distributions
#'   \deqn{\pi_1 \mathrm{Unif}(\mu - a_1, \mu + a_1)
#'   + \ldots + \pi_K \mathrm{Unif}(\mu - a_K, \mu + a_K),}
#'   where parameters \eqn{\pi_k} are estimated and the grid
#'   of (half-)lengths \eqn{(a_1, \ldots, a_K)} is fixed in advance. By
#'   making the grid sufficiently dense, one can obtain an arbitrarily good
#'   approximation. The grid can be specified by the user via parameter
#'   \code{scale}, in which case the argument should be the vector
#'   \eqn{(a_1, \ldots, a_K)}; alternatively, if
#'   \code{scale = "estimate"}, then \code{ebnm} sets the grid via function
#'   \code{\link{ebnm_scale_unimix}}.
#'   Note that \code{ebnm} sets the grid differently from
#'   function \code{\link[ashr]{ash}}. To use the \code{ash} grid, set
#'   \code{scale = "estimate"} and pass in \code{gridmult} as an additional
#'   parameter. See \code{\link[ashr]{ash}} for defaults and details.
#'
#' @export
#'
ebnm_unimodal_symmetric <- function(x,
                                    s = 1,
                                    mode = 0,
                                    scale = "estimate",
                                    g_init = NULL,
                                    fix_g = FALSE,
                                    output = ebnm_output_default(),
                                    control = NULL,
                                    ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "unimodal_symmetric",
                        call = match.call(),
                        ...))
}

#' Solve the EBNM problem using unimodal nonnegative distributions
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   unimodal distributions with support constrained to be greater than the
#'   mode. Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "unimodal_nonnegative"}. For details about the model,
#'   see \code{\link{ebnm}}.
#'
#' @inherit ebnm_unimodal
#'
#' @param scale The nonparametric family of nonnegative unimodal distributions is
#'   approximated via a finite mixture of uniform distributions
#'   \deqn{\pi_1 \mathrm{Unif}(\mu, \mu + a_1)
#'   + \ldots + \pi_K \mathrm{Unif}(\mu, \mu + a_K),}
#'   where parameters \eqn{\pi_k} are estimated and the grid
#'   of lengths \eqn{(a_1, \ldots, a_K)} is fixed in advance. By
#'   making the grid sufficiently dense, one can obtain an arbitrarily good
#'   approximation. The grid can be specified by the user via parameter
#'   \code{scale}, in which case the argument should be the vector of
#'   lengths \eqn{(a_1, \ldots, a_K)}; alternatively, if
#'   \code{scale = "estimate"}, then \code{ebnm} sets the grid via function
#'   \code{\link{ebnm_scale_unimix}}.
#'   Note that \code{ebnm} sets the grid differently from
#'   function \code{\link[ashr]{ash}}. To use the \code{ash} grid, set
#'   \code{scale = "estimate"} and pass in \code{gridmult} as an additional
#'   parameter. See \code{\link[ashr]{ash}} for defaults and details.
#'
#' @export
#'
ebnm_unimodal_nonnegative <- function(x,
                                      s = 1,
                                      mode = 0,
                                      scale = "estimate",
                                      g_init = NULL,
                                      fix_g = FALSE,
                                      output = ebnm_output_default(),
                                      control = NULL,
                                      ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "unimodal_nonnegative",
                        call = match.call(),
                        ...))
}

#' Solve the EBNM problem using unimodal nonpositive distributions
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   unimodal distributions with support constrained to be less than the
#'   mode. Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "unimodal_nonpositive"}. For details about the model,
#'   see \code{\link{ebnm}}.
#'
#' @inherit ebnm_unimodal
#'
#' @param scale The nonparametric family of nonnpositive unimodal distributions is
#'   approximated via a finite mixture of uniform distributions
#'   \deqn{\pi_1 \mathrm{Unif}(\mu - a_1, \mu)
#'   + \ldots + \pi_K \mathrm{Unif}(\mu - a_K, \mu),}
#'   where parameters \eqn{\pi_k} are estimated and the grid
#'   of lengths \eqn{(a_1, \ldots, a_K)} is fixed in advance. By
#'   making the grid sufficiently dense, one can obtain an arbitrarily good
#'   approximation. The grid can be specified by the user via parameter
#'   \code{scale}, in which case the argument should be the vector of
#'   lengths \eqn{(a_1, \ldots, a_K)}; alternatively, if
#'   \code{scale = "estimate"}, then \code{ebnm} sets the grid via function
#'   \code{\link{ebnm_scale_unimix}}.
#'   Note that \code{ebnm} sets the grid differently from
#'   function \code{\link[ashr]{ash}}. To use the \code{ash} grid, set
#'   \code{scale = "estimate"} and pass in \code{gridmult} as an additional
#'   parameter. See \code{\link[ashr]{ash}} for defaults and details.
#'
#' @export
#'
ebnm_unimodal_nonpositive <- function(x,
                                      s = 1,
                                      mode = 0,
                                      scale = "estimate",
                                      g_init = NULL,
                                      fix_g = FALSE,
                                      output = ebnm_output_default(),
                                      control = NULL,
                                      ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "unimodal_nonpositive",
                        call = match.call(),
                        ...))
}

#' Solve the EBNM problem using generalized binary priors
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   nonnegative distributions consisting of mixtures where one component is a
#'   point mass at zero and the other is a truncated normal distribution with
#'   lower bound zero and nonzero mode. Typically, the mode is positive, with
#'   the ratio of the mode to the standard deviation taken to be large, so that
#'   posterior estimates are strongly shrunk towards one of two values (zero or
#'   the mode of the normal component).
#'   Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "generalized_binary"}.
#'   For details, see Liu et al. (2023), cited in \strong{References} below.
#'
#' @inherit ebnm_point_normal
#'
#' @param mode A scalar specifying the mode of the truncated normal component,
#'   or \code{"estimate"} if the mode is to be estimated from the data (the
#'   location of the point mass is fixed at zero).
#'
#' @param scale A scalar specifying the ratio of the (untruncated) standard
#'   deviation of the normal component to its mode. This ratio must be
#'   fixed in advance (i.e., it is not possible to set \code{scale = "estimate"}
#'   when using generalized binary priors).
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link[ashr]{tnormalmix}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{tnormalmix}.
#'
#' @param control A list of control parameters to be passed to function
#'   \code{\link[stats]{optim}}, where \code{method} has been set to
#'   \code{"L-BFGS-B"}.
#'
#' @param ... The following additional arguments act as control parameters for
#'   the outer EM loops in the fitting algorithm. Each loop iteratively updates
#'   parameters \eqn{w} (the
#'   mixture proportion corresponding to the truncated normal component) and
#'   \eqn{\mu} (the mode of the truncated normal component):
#'     \describe{
#'        \item{\code{wlist}}{A vector defining intervals of \eqn{w} for which
#'          optimal solutions will separately be found. For example, if
#'          \code{wlist = c(0, 0.5, 1)}, then two optimal priors will be found:
#'          one such that \eqn{w} is constrained to be less than 0.5 and one
#'          such that it is constrained to be greater than 0.5.}
#'        \item{\code{maxiter}}{A scalar specifying the maximum number of
#'          iterations to perform in each outer EM loop.}
#'        \item{\code{tol}}{A scalar specifying the convergence tolerance
#'          parameter for each outer EM loop.}
#'        \item{\code{mu_init}}{A scalar specifying the initial value of \eqn{\mu}
#'          to be used in each outer EM loop.}
#'        \item{\code{mu_range}}{A vector of length two specifying lower and
#'          upper bounds for possible values of \eqn{\mu}.}
#'      }
#'
#' @references
#' Yusha Liu, Peter Carbonetto, Jason Willwerscheid, Scott A Oakes, Kay F Macleod,
#'   and Matthew Stephens (2023). Dissecting tumor transcriptional heterogeneity
#'   from single-cell RNA-seq data by generalized binary covariance decomposition.
#'   bioRxiv 2023.08.15.553436.
#'
#' @export
#'
ebnm_generalized_binary <- function(x,
                                    s = 1,
                                    mode = "estimate",
                                    scale = 0.1,
                                    g_init = NULL,
                                    fix_g = FALSE,
                                    output = ebnm_output_default(),
                                    control = NULL,
                                    ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "generalized_binary",
                        call = match.call(),
                        ...))
}

#' Solve the EBNM problem using the family of all distributions
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family
#'   of all distributions. Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "npmle"}. For details about the model, see
#'   \code{\link{ebnm}}.
#'
#' @inherit ebnm_normal_scale_mixture
#'
#' @param scale The nonparametric family of all distributions is
#'   approximated via a finite mixture of point masses
#'   \deqn{\pi_1 \delta_{\mu_1} + \ldots + \pi_K \delta_{\mu_K},}
#'   where parameters \eqn{\pi_k} are estimated and the point masses are
#'   evenly spaced over \eqn{(\mu_1, \mu_K)}. By taking a sufficiently dense
#'   grid of point masses, one can obtain an arbitrarily good
#'   approximation. The distance between successive point masses can be
#'   specified by the user via parameter
#'   \code{scale}, in which case the argument should be a scalar specifying the
#'   distance \eqn{d = \mu_2 - \mu_1 = \cdots = \mu_K - \mu_{K - 1}};
#'   alternatively, if \code{scale = "estimate"}, then \code{ebnm} sets the grid
#'   via function \code{\link{ebnm_scale_npmle}}.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. This has
#'   the side effect of fixing the \code{scale} parameter. When supplied,
#'   \code{g_init} should be an object of class \code{\link[ashr]{normalmix}}
#'   or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{normalmix}.
#'
#' @param optmethod Not used by \code{ebnm_npmle}.
#'
#' @export
#'
ebnm_npmle <- function(x,
                       s = 1,
                       scale = "estimate",
                       g_init = NULL,
                       fix_g = FALSE,
                       output = ebnm_output_default(),
                       optmethod = NULL,
                       control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = 0,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = optmethod,
                        control = control,
                        prior_family = "npmle",
                        call = match.call()))
}

#' Solve the EBNM problem using the "deconvolveR" family of distributions
#'
#' Solves the empirical Bayes normal means (EBNM) problem using a non-parametric
#'   exponential family with a natural spline basis.
#'   Like \code{\link{ebnm_npmle}}, there is no unimodal assumption, but whereas
#'   \code{ebnm_npmle} produces spiky estimates for \eqn{g},
#'   \code{ebnm_deconvolver} estimates are much more regular. See
#'   \code{\link[deconvolveR]{deconvolveR-package}} for details and
#'   references. Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "deconvolver"}.
#'
#' @inherit ebnm_npmle
#'
#' @param s Standard errors, which must be uniformly equal to 1 (i.e.,
#'   \code{s = 1}) since the deconvolveR method takes \eqn{z}-scores as input.
#'
#' @param scale A deconvolveR prior is a finite mixture of point masses
#'   \deqn{\pi_1 \delta_{\mu_1} + \ldots + \pi_K \delta_{\mu_K},}
#'   where parameters \eqn{\pi_k} are estimated and the point masses are
#'   evenly spaced over \eqn{(\mu_1, \mu_K)}.The distance between successive
#'   point masses can be specified by the user via parameter
#'   \code{scale}, in which case the argument should be a scalar specifying the
#'   distance \eqn{d = \mu_2 - \mu_1 = \cdots = \mu_K - \mu_{K - 1}};
#'   alternatively, if \code{scale = "estimate"}, then \code{ebnm} sets the grid
#'   via function \code{\link{ebnm_scale_npmle}}.
#'
#' @param control A list of control parameters to be passed to optimization
#'   function \code{\link[stats]{nlm}}.
#'
#' @param ... Additional parameters to be passed to function
#'   \code{\link[deconvolveR]{deconv}} in package \code{deconvolveR}.
#'
#' @export
#'
ebnm_deconvolver <- function(x,
                             s = 1,
                             scale = "estimate",
                             g_init = NULL,
                             fix_g = FALSE,
                             output = ebnm_output_default(),
                             control = NULL,
                             ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = 0,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        control = control,
                        prior_family = "deconvolver",
                        call = match.call(),
                        ...))
}

#' Solve the EBNM problem using a flat prior
#'
#' Solves the empirical Bayes normal means (EBNM) problem using a
#'   "non-informative" improper uniform prior, which yields posteriors
#'   \deqn{\theta_j | x_j, s_j \sim N(x_j, s_j^2).} Identical to function
#'   \code{\link{ebnm}} with argument \code{prior_family = "flat"}. For details
#'   about the model, see \code{\link{ebnm}}.
#'
#' @inherit ebnm_point_normal
#'
#' @param g_init Not used by \code{ebnm_flat}, but included for consistency
#'   with other \code{ebnm} functions.
#'
#' @param fix_g Not used by \code{ebnm_flat}, but included for consistency
#'   with other \code{ebnm} functions.
#'
#' @export
#'
ebnm_flat <- function(x,
                      s = 1,
                      g_init = NULL,
                      fix_g = FALSE,
                      output = ebnm_output_default()) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = 0,
                        scale = 0,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        control = NULL,
                        prior_family = "flat",
                        call = match.call()))
}

#' Solve the EBNM problem using a point mass prior
#'
#' Solves the empirical Bayes normal means (EBNM) problem using the family of
#'   point masses \eqn{\delta_\mu}. Posteriors are simply point masses at \eqn{\mu}.
#'   Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "point_mass"}. For details about the model, see
#'   \code{\link{ebnm}}.
#'
#' @inherit ebnm_point_normal
#'
#' @param mode A scalar specifying the location of the point mass or
#'   \code{"estimate"} if the location is to be estimated from the data.
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. When
#'   supplied, \code{g_init} should be an object of class
#'   \code{\link[ashr]{normalmix}} or an \code{ebnm} object in which the fitted
#'   prior is an object of class \code{normalmix}.
#'
#' @export
#'
ebnm_point_mass <- function(x,
                            s = 1,
                            mode = 0,
                            g_init = NULL,
                            fix_g = FALSE,
                            output = ebnm_output_default()) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = 0,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = NULL,
                        prior_family = "point_mass",
                        call = match.call()))
}

#' Solve the EBNM problem using an ash family of distributions
#'
#' A wrapper to function \code{\link[ashr]{ash}} in package \code{ashr}.
#'   Identical to function \code{\link{ebnm}} with argument
#'   \code{prior_family = "ash"}.
#'
#' @inherit ebnm_normal_scale_mixture
#'
#' @param mode Passed to \code{\link[ashr]{ash}} as parameter \code{mode}.
#'
#' @param scale Passed to \code{\link[ashr]{ash}} as parameter \code{mixsd}.
#'
#' @param g_init Passed to \code{\link[ashr]{ash}} as parameter \code{g}.
#'
#' @param fix_g Passed to \code{\link[ashr]{ash}} as parameter \code{fixg}.
#'
#' @param control Passed to \code{\link[ashr]{ash}} as parameter \code{control}.
#'
#' @param ... Additional parameters to be passed to \code{\link[ashr]{ash}}.
#'
#' @export
#'
ebnm_ash <- function(x,
                     s = 1,
                     mode = 0,
                     scale = "estimate",
                     g_init = NULL,
                     fix_g = FALSE,
                     output = ebnm_output_default(),
                     control = NULL,
                     ...) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = NULL,
                        control = control,
                        prior_family = "ash",
                        call = match.call(),
                        ...))
}
