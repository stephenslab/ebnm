#' Solve the EBNM problem
#'
#' Solves the empirical Bayes normal means problem using a specified family of
#'   priors.
#'
#' Given vectors of data \code{x} and standard errors \code{s},
#'   solve the "empirical Bayes normal means" (EBNM) problem for various
#'   choices of prior family.
#'   The model is \deqn{x_j | \theta_j, s_j \sim N(\theta_j, s_j^2)}
#'   \deqn{\theta_j | s_j \sim g \in G} where the distribution \eqn{g} is to
#'   be estimated.
#'   The distribution \eqn{g} is referred to as the "prior distribution" for
#'   \eqn{\theta}
#'   and \eqn{G} is a specified family of prior distributions. Several options
#'   for \eqn{G} are implemented, some parametric and others non-parametric;
#'   see below for examples.
#'
#'   Solving the EBNM problem involves
#'   two steps. First, estimate \eqn{g \in  G} via maximum marginal likelihood,
#'   yielding an estimate \deqn{\hat{g}:= \arg\max_{g \in G} L(g)} where
#'   \deqn{L(g):= \prod_j \int p(x_j | \theta_j, s_j)  g(d\theta_j)}
#'   Second, compute the posterior distributions
#'   \eqn{p(\theta_j | x_j, s_j, \hat{g})} and/or summaries
#'   such as posterior means and posterior second moments.
#'
#'   The prior families that have been implemented include:
#'     \describe{
#'       \item{\code{point_normal}}{The family of mixtures where one
#'         component is a point mass at \eqn{\mu} and the other is a normal
#'         distribution centered at \eqn{\mu}.}
#'       \item{\code{point_laplace}}{The family of mixtures where one
#'         component is a point mass at zero and the other is a
#'         double-exponential distribution.}
#'       \item{\code{point_exponential}}{The family of mixtures where one
#'         component is a point mass at zero and the other is a
#'         (nonnegative) exponential distribution.}
#'       \item{\code{normal}}{The family of normal distributions.}
#'       \item{\code{horseshoe}}{The family of \link{horseshoe} distributions.}
#'       \item{\code{normal_scale_mixture}}{The family of scale mixtures of
#'         normals.}
#'       \item{\code{unimodal}}{The family of all unimodal distributions.}
#'       \item{\code{unimodal_symmetric}}{The family of symmetric unimodal
#'         distributions.}
#'       \item{\code{unimodal_nonnegative}}{The family of unimodal
#'         distributions with support constrained to be greater than the mode.}
#'       \item{\code{unimodal_nonpositive}}{The family of unimodal
#'         distributions with support constrained to be less than the mode.}
#'       \item{\code{npmle}}{The family of all distributions.}
#'       \item{\code{deconvolver}}{A non-parametric exponential family with
#'         a natural spline basis. Like \code{npmle}, there is no unimodal
#'         assumption, but whereas \code{npmle} produces spiky estimates for
#'         \eqn{g}, \code{deconvolver} estimates are much more regular. See
#'         Narasimhan and Efron (2020) for details.}
#'     }
#'
#' @param x A vector of observations. Missing observations (\code{NA}s) are
#'   not allowed.
#'
#' @param s A vector of standard errors (or a scalar if all are equal).
#'   Standard errors may not be exactly zero, and
#'   missing standard errors are not allowed. Two prior families have
#'   additional restrictions: when \code{prior_family = "horseshoe"}, errors
#'   must be homoskedastic; and since function
#'   \code{\link[deconvolveR]{deconv}} in package \code{deconvolveR} takes
#'   z-scores, it must be true that \code{s = 1} when
#'   \code{prior_family = "deconvolver"}.
#'
#' @param prior_family A character string that specifies the prior family
#'   \eqn{G}. See "Details" below.
#'
#' @param mode A scalar specifying the mode of the prior \eqn{g} or
#'   \code{"estimate"} if the mode is to be estimated from the data. This
#'   parameter is ignored when \code{prior_family} is \code{"npmle"} or
#'   \code{"deconvolver"}.
#'
#' @param scale A scalar or vector specifying the scale parameter(s) of the
#'   prior or \code{"estimate"} if the scale parameters are to be estimated
#'   from the data. The precise interpretation of \code{scale} depends on
#'   \code{prior_family}. For normal and point-normal families, it is a scalar
#'   specifying the standard deviation of the normal component. For
#'   point-Laplace and point-exponential families, it is a scalar specifying
#'   the scale parameter of the Laplace or exponential component. For horseshoe
#'   families, it corresponds to \eqn{s\tau} in the usual parametrization of
#'   the \code{\link{horseshoe}} distribution. For \code{"npmle"} and
#'   \code{"deconvolver"}, it is a scalar specifying the distance between support
#'   points. For all other prior families, which are implemented using the function
#'   \code{\link[ashr]{ash}} in package \code{ashr}, it is a vector specifying
#'   the parameter \code{mixsd} to be passed to \code{ash} or \code{"estimate"}
#'   if \code{mixsd} is to be chosen by \code{ebnm}. (Note that \code{ebnm} chooses
#'   \code{mixsd} differently from \code{ashr}. To use the \code{ashr} grid, set
#'   \code{scale = "estimate"} and pass in \code{gridmult}, \code{pointmass}, or
#'   \code{method} as an additional \code{ash} parameter. See
#'   \code{\link[ashr]{ash}} for defaults and details.)
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. For
#'   non-parametric priors, this has the side effect of fixing the \code{mode}
#'   and \code{scale} parameters. If \code{g_init} is supplied, it should be
#'   an object of class \code{\link[ashr]{normalmix}} for prior families
#'   \code{normal}, \code{point_normal}, \code{normal_scale_mixture},
#'   \code{npmle}, and \code{deconvolver}; class \code{\link{laplacemix}} for
#'   point-Laplace families; class \code{\link{gammamix}} for point-exponential
#'   families; class \code{\link{horseshoe}} for horseshoe families; or class
#'   \code{\link[ashr]{unimix}} for \code{unimodal_} families.
#'
#' @param fix_g If \code{TRUE}, fix the prior \eqn{g} at \code{g_init} instead
#'   of estimating it.
#'
#' @param output A character vector indicating which values are to be returned.
#'   Function \code{output_default()} provides the default return values, while
#'   \code{output_all()} lists all possible return values. See "Value" below.
#'
#' @param optmethod A string specifying which optimization function is to be
#'   used. Options include \code{"nlm"}, \code{"lbfgsb"} (which calls
#'   \code{optim} with \code{method = "L-BFGS-B"}), and \code{"trust"} (which
#'   calls into package \code{trust}). Other options are \code{"nohess_nlm"},
#'   \code{"nograd_nlm"}, and \code{"nograd_lbfgsb"}, which use numerical
#'   approximations rather than exact expressions for the Hessian and (for
#'   the latter two) the gradient. The default option is \code{"nohess_nlm"}.
#'   Since all non-parametric families rely upon external packages, this parameter
#'   is only available for parametric families (point-normal, point-Laplace,
#'   point-exponential, and normal).
#'
#' @param control A list of control parameters to be passed to the optimization
#'   function. \code{\link[stats]{optimize}} is used for
#'   \code{prior_family = "normal"} and \code{prior_family = "horseshoe"},
#'   while \code{\link[stats]{nlm}} is used for parametric families unless
#'   parameter \code{optmethod} specifies otherwise. \code{\link[stats]{nlm}} is
#'   also used for \code{prior_family = "deconvolver"}.
#'   For ash families (including \code{normal_scale_mixture}, \code{npmle}, and
#'   all \code{unimodal_} families), function \code{\link[mixsqp]{mixsqp}} in
#'   package \code{mixsqp} is the default.
#'
#' @param ... Additional parameters. When \code{prior_family = "ash"} or when
#'   a \code{unimodal_} prior is used, these parameters are passed to
#'   function \code{\link[ashr]{ash}} in package \code{ashr}. When
#'   \code{prior_family = "deconvolver"}, they are passed to function
#'   \code{\link[deconvolveR]{deconv}} in package \code{deconvolveR}. Although it
#'   does not call into \code{ashr}, \code{normal_scale_mixture} accepts
#'   parameter \code{gridmult} for purposes of comparison. When \code{gridmult}
#'   is set, an \code{ashr}-style grid will be used instead of the default
#'   \code{ebnm} grid. In all other cases, additional parameters are ignored.
#'
#' @return An \code{ebnm} object. Depending on the argument to \code{output}, the
#'   object is a list containing elements:
#'     \describe{
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, standard deviations, and second moments; local false sign
#'         rates).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}} (an object of
#'         class \code{\link[ashr]{normalmix}}, \code{\link{laplacemix}},
#'         \code{\link{gammamix}}, \code{\link[ashr]{unimix}}, or
#'         \code{\link{horseshoe}}).}
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained
#'         \eqn{L(\hat{g})}.}
#'       \item{\code{posterior_sampler}}{A function that can be used to
#'         produce samples from the posterior. For all prior families other
#'         than the horseshoe, it takes a single parameter
#'         \code{nsamp}, the number of posterior samples to return per
#'         observation. Since \code{ebnm_horseshoe} returns an MCMC sampler,
#'         it additionally takes parameter \code{burn}, the number of burn-in
#'         samples to discard.}
#'      }
#'
#' @seealso Calling functions \code{\link{ebnm_point_normal}},
#'   \code{\link{ebnm_point_laplace}},
#'   \code{\link{ebnm_point_exponential}}, \code{\link{ebnm_normal}},
#'   \code{\link{ebnm_horseshoe}},
#'   \code{\link{ebnm_normal_scale_mixture}}, \code{\link{ebnm_unimodal}},
#'   \code{\link{ebnm_unimodal_symmetric}},
#'   \code{\link{ebnm_unimodal_nonnegative}},
#'   \code{\link{ebnm_unimodal_nonpositive}}, \code{\link{ebnm_npmle}},
#'   \code{\link{ebnm_deconvolver}}, and \code{\link{ebnm_ash}}
#'   is equivalent to calling \code{ebnm} with \code{prior_family} set
#'   accordingly.
#'
#' @examples
#' theta <- c(rep(0, 100), rexp(100))
#' s <- 1
#' x <- theta + rnorm(200, 0, s)
#'
#' # The following are equivalent:
#' pn.res <- ebnm(x, s, prior_family = "point_normal")
#' pn.res <- ebnm_point_normal(x, s)
#'
#' # Inspect results:
#' pn.res$log_likelihood
#' plot(x, pn.res$posterior$mean)
#'
#' # Fix the scale parameter:
#' pl.res <- ebnm_point_laplace(x, s, scale = 1)
#' pl.res$fitted_g$scale
#'
#' # Estimate the mode:
#' normal.res <- ebnm_normal(x, s, mode = "estimate")
#' normal.res$fitted_g$mean
#'
#' # Use an initial g (this fixes mode and scale for ash priors):
#' normalmix.res <- ebnm_normal_scale_mixture(x, s, g_init = pn.res$fitted_g)
#'
#' # Fix g and get more output:
#' g_init <- pn.res$fitted_g
#' pn.res <- ebnm_point_normal(x, s, g_init = g_init, fix_g = TRUE,
#'                             output = "posterior_sampler")
#' pn.res <- ebnm_point_normal(x, s, g_init = g_init, fix_g = TRUE,
#'                             output = output_all())
#'
#' # Examples of usage of control parameter:
#' #  point_normal uses nlm:
#' pn.res <- ebnm_point_normal(x, s, control = list(print.level = 1))
#' #  unimodal uses mixsqp:
#' unimodal.res <- ebnm_unimodal(x, s, control = list(verbose = TRUE))
#'
#' @export
#'
ebnm <- function(x,
                 s = 1,
                 prior_family = c("point_normal",
                                  "point_laplace",
                                  "point_exponential",
                                  "normal",
                                  "horseshoe",
                                  "normal_scale_mixture",
                                  "unimodal",
                                  "unimodal_symmetric",
                                  "unimodal_nonnegative",
                                  "unimodal_nonpositive",
                                  "npmle",
                                  "deconvolver",
                                  "ash"),
                 mode = 0,
                 scale = "estimate",
                 g_init = NULL,
                 fix_g = FALSE,
                 output = output_default(),
                 optmethod = NULL,
                 control = NULL,
                 ...) {
  prior_family <- match.arg(prior_family)

  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        optmethod = optmethod,
                        control = control,
                        prior_family = prior_family,
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
                           optmethod,
                           control,
                           prior_family,
                           call,
                           ...) {
  check_args(x, s, g_init, fix_g, output, mode)
  mode <- handle_mode_parameter(mode)
  scale <- handle_scale_parameter(scale)
  if (is.null(control)) {
    control <- list()
  }
  if (lfsr_arg_str() %in% output && mode != 0) {
    warning("Since they're not well defined for nonzero modes, local false ",
            "sign rates won't be returned.")
    output <- setdiff(output, lfsr_arg_str())
  }

  if (prior_family == "point_normal") {
    retlist <- parametric_workhorse(x = x,
                                    s = s,
                                    mode = mode,
                                    scale = scale,
                                    pointmass = TRUE,
                                    g_init = g_init,
                                    fix_g = fix_g,
                                    output = output,
                                    optmethod = optmethod,
                                    control = control,
                                    checkg_fn = pn_checkg,
                                    initpar_fn = pn_initpar,
                                    scalepar_fn = pn_scalepar,
                                    precomp_fn = pn_precomp,
                                    nllik_fn = pn_nllik,
                                    postcomp_fn = pn_postcomp,
                                    summres_fn = pn_summres,
                                    partog_fn = pn_partog,
                                    postsamp_fn = pn_postsamp,
                                    call = call)
  } else if (prior_family == "point_laplace") {
    retlist <- parametric_workhorse(x = x,
                                    s = s,
                                    mode = mode,
                                    scale = scale,
                                    pointmass = TRUE,
                                    g_init = g_init,
                                    fix_g = fix_g,
                                    output = output,
                                    optmethod = optmethod,
                                    control = control,
                                    checkg_fn = pl_checkg,
                                    initpar_fn = pl_initpar,
                                    scalepar_fn = pl_scalepar,
                                    precomp_fn = pl_precomp,
                                    nllik_fn = pl_nllik,
                                    postcomp_fn = pl_postcomp,
                                    summres_fn = pl_summres,
                                    partog_fn = pl_partog,
                                    postsamp_fn = pl_postsamp,
                                    call = call)
  } else if (prior_family == "point_exponential") {
    retlist <- parametric_workhorse(x = x,
                                    s = s,
                                    mode = mode,
                                    scale = scale,
                                    pointmass = TRUE,
                                    g_init = g_init,
                                    fix_g = fix_g,
                                    output = output,
                                    optmethod = optmethod,
                                    control = control,
                                    checkg_fn = pe_checkg,
                                    initpar_fn = pe_initpar,
                                    scalepar_fn = pe_scalepar,
                                    precomp_fn = pe_precomp,
                                    nllik_fn = pe_nllik,
                                    postcomp_fn = pe_postcomp,
                                    summres_fn = pe_summres,
                                    partog_fn = pe_partog,
                                    postsamp_fn = pe_postsamp,
                                    call = call)
  } else if (prior_family == "normal") {
    retlist <- parametric_workhorse(x = x,
                                    s = s,
                                    mode = mode,
                                    scale = scale,
                                    pointmass = FALSE,
                                    g_init = g_init,
                                    fix_g = fix_g,
                                    output = output,
                                    optmethod = optmethod,
                                    control = control,
                                    checkg_fn = pn_checkg,
                                    initpar_fn = pn_initpar,
                                    scalepar_fn = pn_scalepar,
                                    precomp_fn = pn_precomp,
                                    nllik_fn = pn_nllik,
                                    postcomp_fn = pn_postcomp,
                                    summres_fn = pn_summres,
                                    partog_fn = pn_partog,
                                    postsamp_fn = pn_postsamp,
                                    call = call)
  } else if (prior_family == "horseshoe") {
    retlist <- horseshoe_workhorse(x = x,
                                   s = s,
                                   mode = mode,
                                   scale = scale,
                                   g_init = g_init,
                                   fix_g = fix_g,
                                   output = output,
                                   control = control,
                                   call = call)
  } else if (prior_family == "normal_scale_mixture") {
    retlist <- ebnm_normal_mix_workhorse(x = x,
                                         s = s,
                                         mode = mode,
                                         scale = scale,
                                         g_init = g_init,
                                         fix_g = fix_g,
                                         output = output,
                                         control = control,
                                         call = call,
                                         ...)
  } else if (prior_family == "unimodal") {
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
  } else if (prior_family == "unimodal_symmetric") {
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
  } else if (prior_family == "unimodal_nonnegative") {
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
  } else if (prior_family == "unimodal_nonpositive") {
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
  } else if (prior_family == "ash") {
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
  } else if (prior_family == "npmle"
             && !is.null(optmethod) && optmethod == "REBayes") {
    retlist <- rebayes_workhorse(x = x,
                                 s = s,
                                 mode = mode,
                                 scale = scale,
                                 g_init = g_init,
                                 fix_g = fix_g,
                                 output = output,
                                 control = control,
                                 call = call,
                                 ...)
  } else if (prior_family == "npmle") {
    if (is.null(g_init)) {
      if (!is.null(call$mode)) {
        warning("mode parameter is ignored by ebnm_npmle.")
        call$mode <- NULL
      }

      g_init <- init_g_for_npmle(x, s, scale)
      call$scale <- NULL
    }

    retlist <- ebnm_normal_mix_workhorse(x = x,
                                         s = s,
                                         mode = mode,
                                         scale = scale,
                                         g_init = g_init,
                                         fix_g = fix_g,
                                         output = output,
                                         control = control,
                                         call = call,
                                         ...)
  } else if (prior_family == "deconvolver") {
    retlist <- deconvolver_workhorse(x = x,
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

  return(as_ebnm(retlist))
}

check_args <- function(x, s, g_init, fix_g, output, mode) {
  if (!(length(s) %in% c(1, length(x)))) {
    stop("Argument 's' must have either length 1 or the same length as ",
         "argument 'x'.")
  }

  # if (any(is.na(x)) && !all(is.infinite(s[is.na(x)]))) {
  #   stop("All missing observations must have infinite SEs.")
  # }

  if (any(is.na(x))) {
    stop("Missing observations are not allowed.")
  }

  if (any(is.na(s))) {
    stop("Missing standard errors are not allowed.")
  }

  # Remove this check when handling of zero SEs has been implemented for
  #   point-Laplace and normal-mixture priors and issue #84 in ashr has been
  #   fixed.
  if (any(s <= 0)) {
    stop("Standard errors must be positive (and nonzero).")
  }

  if (any(is.infinite(s))) {
    stop("Standard errors cannot be infinite.")
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
