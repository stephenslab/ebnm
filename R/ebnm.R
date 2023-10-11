#' Solve the EBNM problem
#'
#' Solves the empirical Bayes normal means (EBNM) problem using a specified
#'   family of priors. For an article-length introduction to the package, see
#'   Willwerscheid and Stephens (2021), cited in \strong{References} below.
#'
#' Given vectors of data \code{x} and standard errors \code{s}, \code{ebnm}
#'   solves the "empirical Bayes normal means" (EBNM) problem for various
#'   choices of prior family.
#'   The model is \deqn{x_j | \theta_j, s_j \sim N(\theta_j, s_j^2)}
#'   \deqn{\theta_j \sim g \in G,} where \eqn{g}, which is referred to as the
#'   "prior distribution" for \eqn{\theta}, is to be estimated from among
#'   some specified family of prior distributions \eqn{G}. Several options
#'   for \eqn{G} are implemented, some parametric and others non-parametric;
#'   see below for examples.
#'
#'   Solving the EBNM problem involves
#'   two steps. First, \eqn{g \in  G} is estimated via maximum marginal likelihood:
#'   \deqn{\hat{g} := \arg\max_{g \in G} L(g),} where
#'   \deqn{L(g) := \prod_j \int p(x_j | \theta_j, s_j)  g(d\theta_j).}
#'   Second, posterior distributions
#'   \eqn{p(\theta_j | x_j, s_j, \hat{g})} and/or summaries
#'   such as posterior means and posterior second moments are computed.
#'
#'   Implemented prior families include:
#'     \describe{
#'       \item{\code{point_normal}}{The family of mixtures where one
#'         component is a point mass at \eqn{\mu} and the other is a normal
#'         distribution centered at \eqn{\mu}.}
#'       \item{\code{point_laplace}}{The family of mixtures where one
#'         component is a point mass at \eqn{\mu} and the other is a
#'         double-exponential distribution centered at \eqn{\mu}.}
#'       \item{\code{point_exponential}}{The family of mixtures where one
#'         component is a point mass at \eqn{\mu} and the other is a
#'         (nonnegative) exponential distribution with mode \eqn{\mu}.}
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
#'       \item{\code{generalized_binary}}{The family of mixtures where one
#'         component is a point mass at zero and the other is a truncated
#'         normal distribution with lower bound zero and nonzero mode. See
#'         Liu et al. (2023), cited in \strong{References} below.}
#'       \item{\code{npmle}}{The family of all distributions.}
#'       \item{\code{deconvolver}}{A non-parametric exponential family with
#'         a natural spline basis. Like \code{npmle}, there is no unimodal
#'         assumption, but whereas \code{npmle} produces spiky estimates for
#'         \eqn{g}, \code{deconvolver} estimates are much more regular. See
#'         \code{\link[deconvolveR]{deconvolveR-package}} for details and
#'         references.}
#'       \item{\code{flat}}{The "non-informative" improper uniform prior, which
#'         yields posteriors \deqn{\theta_j | x_j, s_j \sim N(x_j, s_j^2).}}
#'       \item{\code{point_mass}}{The family of point masses \eqn{\delta_\mu}.
#'         Posteriors are point masses at \eqn{\mu}.}
#'       \item{\code{ash}}{Calls into function \code{\link[ashr]{ash}} in
#'         package \code{ashr}. Can be used to make direct comparisons of
#'         \code{ebnm} and \code{ashr} implementations of prior families such as
#'         scale mixtures of normals and the NPMLE.}
#'     }
#'
#' @param x A vector of observations. Missing observations (\code{NA}s) are
#'   not allowed.
#'
#' @param s A vector of standard errors (or a scalar if all are equal).
#'   Standard errors may not be exactly zero, and
#'   missing standard errors are not allowed. Two prior families have
#'   additional restrictions: when horseshoe priors are used, errors
#'   must be homoskedastic; and since function
#'   \code{\link[deconvolveR]{deconv}} in package \code{deconvolveR} takes
#'   \eqn{z}-scores, the "deconvolver" family requires that all standard errors
#'   be equal to 1.
#'
#' @param prior_family A character string that specifies the prior family
#'   \eqn{G}. See \strong{Details} below.
#'
#' @param mode A scalar specifying the mode of the prior \eqn{g} or
#'   \code{"estimate"} if the mode is to be estimated from the data. This
#'   parameter is ignored by the NPMLE, the \code{deconvolveR} family,
#'   and the improper uniform (or "flat") prior. For generalized binary priors,
#'   which are bimodal, the mode parameter specifies the mode of the truncated
#'   normal component (the location of the point mass is fixed at zero).
#'
#' @param scale A scalar or vector specifying the scale parameter(s) of the
#'   prior or \code{"estimate"} if the scale parameters are to be estimated
#'   from the data. This parameter is ignored by the flat prior and the family
#'   of point mass priors.
#'
#'   The interpretation of \code{scale} depends on the prior
#'   family. For normal and point-normal families, it is a scalar
#'   specifying the standard deviation of the normal component. For
#'   point-Laplace and point-exponential families, it is a scalar specifying
#'   the scale parameter of the Laplace or exponential component. For the horseshoe
#'   family, it corresponds to \eqn{s\tau} in the usual parametrization of
#'   the \code{\link{horseshoe}} distribution. For the family of generalized
#'   binary priors, it specifies the ratio of the (untruncated) standard
#'   deviation of the normal component to its mode.
#'   This ratio must be fixed in advance (i.e., argument \code{"estimate"} is
#'   unavailable for generalized binary priors). For the NPMLE and \code{deconvolveR}
#'   prior family, \code{scale} is a scalar specifying the distance between
#'   successive means in the grid of point masses or normal distributions
#'   used to estimate \eqn{g}.
#'   For all other prior families, which are implemented using the function
#'   \code{\link[ashr]{ash}} in package \code{ashr}, it is a vector specifying
#'   the parameter \code{mixsd} to be passed to \code{ash} or \code{"estimate"}
#'   if \code{mixsd} is to be chosen by \code{ebnm}. (Note that \code{ebnm} chooses
#'   \code{mixsd} differently from \code{ash}: see functions
#'   \code{\link{ebnm_scale_normalmix}}, \code{\link{ebnm_scale_unimix}}, and
#'   \code{\link{ebnm_scale_npmle}} for details. To use the \code{ash} grid, set
#'   \code{scale = "estimate"} and pass in \code{gridmult} as an additional
#'   parameter. See \code{\link[ashr]{ash}} for defaults and details.)
#'
#' @param g_init The prior distribution \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. For
#'   non-parametric priors, this has the side effect of fixing the \code{mode}
#'   and \code{scale} parameters. If \code{g_init} is supplied, it should be
#'   an object of class \code{\link[ashr]{normalmix}} for normal, point-normal,
#'   scale mixture of normals, and \code{deconvolveR} prior families, as well as
#'   for the NPMLE; class \code{\link{laplacemix}} for
#'   point-Laplace families; class \code{\link{gammamix}} for point-exponential
#'   families; class \code{\link{horseshoe}} for horseshoe families; class
#'   \code{\link[ashr]{unimix}} for \code{unimodal_} families; or class
#'   \code{\link[ashr]{tnormalmix}} for generalized binary priors. An object of
#'   class \code{ebnm} can also be supplied as argument, provided that field
#'   \code{fitted_g} contains a prior of the correct class (see
#'   \strong{Examples} below).
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
#'   used. Since all non-parametric families rely upon external packages, this
#'   parameter is only available for parametric families (point-normal,
#'   point-Laplace, point-exponential, and normal). Options include \code{"nlm"},
#'   \code{"lbfgsb"} (which calls
#'   \code{optim} with \code{method = "L-BFGS-B"}), and \code{"trust"} (which
#'   calls into package \code{trust}). Other options are \code{"nohess_nlm"},
#'   \code{"nograd_nlm"}, and \code{"nograd_lbfgsb"}, which use numerical
#'   approximations rather than exact expressions for the Hessian and (for
#'   the latter two) the gradient. The default option is \code{"nohess_nlm"}.
#'
#'
#' @param control A list of control parameters to be passed to the optimization
#'   function. \code{\link[stats]{optimize}} is used for normal and horseshoe
#'   prior families, while \code{\link[stats]{nlm}} is used for parametric
#'   families unless parameter \code{optmethod} specifies otherwise.
#'   \code{\link[stats]{nlm}} is also used for the \code{deconvolveR} prior family.
#'   For ash families (including scale mixtures of normals, the NPMLE, and
#'   all \code{unimodal_} families), function \code{\link[mixsqp]{mixsqp}} in
#'   package \code{mixsqp} is the default. For generalized binary priors,
#'   function \code{\link[stats]{optim}} is used with \code{method = "L-BFGS-B"}.
#'
#' @param ... Additional parameters. When a \code{unimodal_} prior family is used,
#'   these parameters are passed to function \code{\link[ashr]{ash}} in package
#'   \code{ashr}. Although it
#'   does not call into \code{ashr}, the scale mixture of normals family accepts
#'   parameter \code{gridmult} for purposes of comparison. When \code{gridmult}
#'   is set, an \code{ashr}-style grid will be used instead of the default
#'   \code{ebnm} grid. When the "deconvolver" family is used, additional
#'   parameters are passed to function \code{\link[deconvolveR]{deconv}} in
#'   package \code{deconvolveR}. Families of generalized binary priors take several
#'   additional parameters; see \code{\link{ebnm_generalized_binary}}. In all
#'   other cases, additional parameters are ignored.
#'
#' @return An \code{ebnm} object. Depending on the argument to \code{output}, the
#'   object is a list containing elements:
#'     \describe{
#'       \item{\code{data}}{A data frame containing the observations \code{x}
#'         and standard errors \code{s}.}
#'       \item{\code{posterior}}{A data frame of summary results (posterior
#'         means, standard deviations, second moments, and local false sign
#'         rates).}
#'       \item{\code{fitted_g}}{The fitted prior \eqn{\hat{g}} (an object of
#'         class \code{\link[ashr]{normalmix}}, \code{\link{laplacemix}},
#'         \code{\link{gammamix}}, \code{\link[ashr]{unimix}},
#'         \code{\link[ashr]{tnormalmix}}, or \code{\link{horseshoe}}).}
#'       \item{\code{log_likelihood}}{The optimal log likelihood attained,
#'         \eqn{L(\hat{g})}.}
#'       \item{\code{posterior_sampler}}{A function that can be used to
#'         produce samples from the posterior. For all prior families other
#'         than the horseshoe, the sampler takes a single parameter
#'         \code{nsamp}, the number of posterior samples to return per
#'         observation. Since \code{ebnm_horseshoe} returns an MCMC sampler,
#'         it additionally takes parameter \code{burn}, the number of burn-in
#'         samples to discard.}
#'      }
#'    S3 methods \code{coef}, \code{confint}, \code{fitted}, \code{logLik},
#'    \code{nobs}, \code{plot}, \code{predict}, \code{print}, \code{quantile},
#'    \code{residuals}, \code{simulate}, \code{summary}, and \code{vcov}
#'    have been implemented for \code{ebnm} objects. For details, see the
#'    respective help pages, linked below under \strong{See Also}.
#'
#' @references
#' Jason Willwerscheid and Matthew Stephens (2021).
#'   \code{ebnm}: an \code{R} Package for solving the empirical Bayes
#'   normal means problem using a variety of prior families. arXiv:2110.00152.
#'
#' Yusha Liu, Peter Carbonetto, Jason Willwerscheid, Scott A Oakes, Kay F Macleod,
#'   and Matthew Stephens (2023). Dissecting tumor transcriptional heterogeneity
#'   from single-cell RNA-seq data by generalized binary covariance decomposition.
#'   bioRxiv 2023.08.15.553436.

#'
#' @seealso A plotting method is available for \code{ebnm} objects: see
#'   \code{\link{plot.ebnm}}.
#'
#'   For other methods, see \code{\link{coef.ebnm}}, \code{\link{confint.ebnm}},
#'    \code{\link{fitted.ebnm}}, \code{\link{logLik.ebnm}},
#'    \code{\link{nobs.ebnm}}, \code{\link{predict.ebnm}},
#'    \code{\link{print.ebnm}}, \code{\link{print.summary.ebnm}},
#'    \code{\link{quantile.ebnm}}, \code{\link{residuals.ebnm}},
#'    \code{\link{simulate.ebnm}},
#'    \code{\link{summary.ebnm}}, and \code{\link{vcov.ebnm}}.
#'
#'   Calling into functions \code{\link{ebnm_point_normal}},
#'   \code{\link{ebnm_point_laplace}},
#'   \code{\link{ebnm_point_exponential}}, \code{\link{ebnm_normal}},
#'   \code{\link{ebnm_horseshoe}},
#'   \code{\link{ebnm_normal_scale_mixture}}, \code{\link{ebnm_unimodal}},
#'   \code{\link{ebnm_unimodal_symmetric}},
#'   \code{\link{ebnm_unimodal_nonnegative}},
#'   \code{\link{ebnm_unimodal_nonpositive}},
#'   \code{\link{ebnm_generalized_binary}}, \code{\link{ebnm_npmle}},
#'   \code{\link{ebnm_deconvolver}}, \code{\link{ebnm_flat}},
#'   \code{\link{ebnm_point_mass}}, and \code{\link{ebnm_ash}}
#'   is equivalent to calling into \code{ebnm} with \code{prior_family} set
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
#' logLik(pn.res)
#' plot(pn.res)
#'
#' # Fix the scale parameter:
#' pl.res <- ebnm_point_laplace(x, s, scale = 1)
#'
#' # Estimate the mode:
#' normal.res <- ebnm_normal(x, s, mode = "estimate")
#' plot(normal.res) # posterior means shrink to a value different from zero
#'
#' # Use an initial g (this fixes mode and scale for ash priors):
#' normalmix.res <- ebnm_normal_scale_mixture(x, s, g_init = pn.res)
#'
#' # Fix g and get different output (including a posterior sampler):
#' pn.res <- ebnm_point_normal(x, s, g_init = pn.res, fix_g = TRUE,
#'                             output = ebnm_output_all())
#'
#' # Sample from the posterior:
#' pn.samp <- simulate(pn.res, nsim = 100)
#'
#' # Quantiles and HPD confidence intervals can be obtained via sampling:
#' set.seed(1)
#' pn.quantiles <- quantile(pn.res, probs = c(0.1, 0.9))
#' pn.quantiles[1:5, ]
#' confint(pn.res, level = 0.8, parm = 1:5)
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
                                  "generalized_binary",
                                  "npmle",
                                  "deconvolver",
                                  "flat",
                                  "point_mass",
                                  "ash"),
                 mode = 0,
                 scale = "estimate",
                 g_init = NULL,
                 fix_g = FALSE,
                 output = ebnm_output_default(),
                 optmethod = NULL,
                 control = NULL,
                 ...) {

  # Allow the prior family to remain unspecified if it can be inferred from
  #   g_init. This functionality is not currently documented, and has thus not
  #   been tested exhaustively, but it is needed by method predict.
  if (missing(prior_family) && !is.null(g_init)) {
    prior_family <- infer_prior_family(g_init)
  }

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
  if (inherits(g_init, "ebnm")) {
    g_init <- g_init[[g_ret_str()]]
  }

  check_args(x, s, g_init, fix_g, output, mode)
  s <- handle_standard_errors(x, s)

  # Convenience function:
  if (prior_family == "point_mass") {
    prior_family <- "normal"
    if (is.null(g_init) && !is.null(call$scale)) {
      warning("scale parameter is ignored by ebnm_point_mass.")
      call$scale <- NULL
    }
    scale <- 0

    if (!is.null(g_init) &&
        (!inherits(g_init, "normalmix") ||
          !(length(g_init$pi) == 1) ||
          !(g_init$sd == 0))) {
      stop("If specified, g_init must be an object of class normalmix consisting ",
           "of a single mixture component with standard deviation zero.")
    }
  }

  mode <- handle_mode_parameter(mode)
  scale <- handle_scale_parameter(scale)

  if (is.null(control)) {
    control <- list()
  }
  if (lfsr_arg_str() %in% output
      && prior_family != "generalized_binary" # gb mode is always zero
      && mode != 0) {
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
                                    call = call,
                                    ...)
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
                                    call = call,
                                    ...)
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
                                    call = call,
                                    ...)
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
                                    call = call,
                                    ...)
  } else if (prior_family == "horseshoe") {
    retlist <- horseshoe_workhorse(x = x,
                                   s = s,
                                   mode = mode,
                                   scale = scale,
                                   g_init = g_init,
                                   fix_g = fix_g,
                                   output = output,
                                   control = control,
                                   call = call,
                                   ...)
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
    # This option is not documented but is used in the paper.
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
  } else if (prior_family == "generalized_binary") {
    # Generalized binary priors have different defaults from other families:
    if (is.null(call$mode)) {
      mode <- "estimate"
    }
    if (is.null(call$scale)) {
      scale <- 0.1
    }
    retlist <- gb_workhorse(x = x,
                            s = s,
                            mode = mode,
                            scale = scale,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            control = control,
                            call = call,
                            ...)
  } else if (prior_family == "flat") {
    if (!(is.null(call$mode) && is.null(call$scale))) {
      warning("mode and scale parameters are ignored by ebnm_flat.")
      call$mode <- NULL
      call$scale <- NULL
    }
    if (!(is.null(call$g_init) && is.null(call$fixg))) {
      warning("g_init and fixg parameters are ignored by ebnm_flat.")
      call$g_init <- NULL
      call$fixg <- NULL
    }
    retlist <- flat_workhorse(x = x,
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

  return(as_ebnm(retlist, call))
}

check_args <- function(x, s, g_init, fix_g, output, mode) {
  if (!(length(s) %in% c(1, length(x)))) {
    stop("Argument 's' must have either length 1 or the same length as ",
         "argument 'x'.")
  }

  if (any(is.na(x))) {
    stop("Missing observations are not allowed.")
  }

  # if (any(is.na(x)) && !all(is.infinite(s[is.na(x)]))) {
  #   stop("All missing observations must have infinite SEs.")
  # }

  if (any(is.na(s))) {
    stop("Missing standard errors are not allowed.")
  }

  if (any(is.infinite(s))) {
    stop("Standard errors cannot be infinite.")
  }

  if (fix_g && is.null(g_init)) {
    stop("If g is fixed, then an initial g must be provided.")
  }

  if (!all(output %in% ebnm_output_all())) {
    stop("Invalid argument to output. See function ebnm_output_all() for a list ",
         "of valid outputs.")
  }
}

handle_standard_errors <- function(x, s) {
  # Remove this check when handling of zero SEs has been implemented for
  #   point-Laplace and normal-mixture priors and issue #84 in ashr has been
  #   fixed.
  if (any(s <= 0)) {
    if (any(s > 0)) {
      min.s <- min(s[s > 0])
    } else {
      min.s <- Inf
    }
    s[s <= 0] <- min(min.s, (max(x) - min(x)) * sqrt(.Machine$double.eps))
    warning("Nonpositive SEs have been replaced by small positive SEs.")
  }

  return(s)
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
