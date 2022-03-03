#' Solve the EBNM problem using a specified family of priors
#'
#' Each of the functions listed below solves the empirical Bayes normal means
#'   (EBNM) problem using a specified family of priors. Calling function
#'   \code{ebnm_xxx} is equivalent to calling function \code{ebnm} with argument
#'   \code{prior_family = "xxx"}. For details about the model, see
#'   \code{\link{ebnm}} or the paper cited in \strong{References} below.
#'
#'   Implemented prior families include:
#'     \describe{
#'       \item{\code{ebnm_point_normal}}{The family of mixtures where one
#'         component is a point mass at \eqn{\mu} and the other is a normal
#'         distribution centered at \eqn{\mu}.}
#'       \item{\code{ebnm_point_laplace}}{The family of mixtures where one
#'         component is a point mass at zero and the other is a
#'         double-exponential distribution.}
#'       \item{\code{ebnm_point_exponential}}{The family of mixtures where one
#'         component is a point mass at zero and the other is a
#'         (nonnegative) exponential distribution.}
#'       \item{\code{ebnm_normal}}{The family of normal distributions.}
#'       \item{\code{ebnm_horseshoe}}{The family of \link{horseshoe} distributions.}
#'       \item{\code{ebnm_normal_scale_mixture}}{The family of scale mixtures of
#'         normals.}
#'       \item{\code{ebnm_unimodal}}{The family of all unimodal distributions.}
#'       \item{\code{ebnm_unimodal_symmetric}}{The family of symmetric unimodal
#'         distributions.}
#'       \item{\code{ebnm_unimodal_nonnegative}}{The family of unimodal
#'         distributions with support constrained to be greater than the mode.}
#'       \item{\code{ebnm_unimodal_nonpositive}}{The family of unimodal
#'         distributions with support constrained to be less than the mode.}
#'       \item{\code{ebnm_npmle}}{The family of all distributions.}
#'       \item{\code{ebnm_deconvolver}}{A non-parametric exponential family with
#'         a natural spline basis. Like \code{npmle}, there is no unimodal
#'         assumption, but whereas \code{npmle} produces spiky estimates for
#'         \eqn{g}, \code{deconvolver} estimates are much more regular. See
#'         \code{\link[deconvolveR]{deconvolveR-package}} for details and
#'         references.}
#'     }
#'
#' @inherit ebnm
#'
#' @seealso \code{\link{ebnm}}
#'
#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_point_normal <- function(x,
                              s = 1,
                              mode = 0,
                              scale = "estimate",
                              g_init = NULL,
                              fix_g = FALSE,
                              output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_point_laplace <- function(x,
                               s = 1,
                               mode = 0,
                               scale = "estimate",
                               g_init = NULL,
                               fix_g = FALSE,
                               output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_point_exponential <- function(x,
                                   s = 1,
                                   mode = 0,
                                   scale = "estimate",
                                   g_init = NULL,
                                   fix_g = FALSE,
                                   output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_normal <- function(x,
                        s = 1,
                        mode = 0,
                        scale = "estimate",
                        g_init = NULL,
                        fix_g = FALSE,
                        output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_horseshoe <- function(x,
                           s = 1,
                           scale = "estimate",
                           g_init = NULL,
                           fix_g = FALSE,
                           output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_normal_scale_mixture <- function(x,
                                      s = 1,
                                      mode = 0,
                                      scale = "estimate",
                                      g_init = NULL,
                                      fix_g = FALSE,
                                      output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_unimodal <- function(x,
                          s = 1,
                          mode = 0,
                          scale = "estimate",
                          g_init = NULL,
                          fix_g = FALSE,
                          output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_unimodal_symmetric <- function(x,
                                    s = 1,
                                    mode = 0,
                                    scale = "estimate",
                                    g_init = NULL,
                                    fix_g = FALSE,
                                    output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_unimodal_nonnegative <- function(x,
                                      s = 1,
                                      mode = 0,
                                      scale = "estimate",
                                      g_init = NULL,
                                      fix_g = FALSE,
                                      output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_unimodal_nonpositive <- function(x,
                                      s = 1,
                                      mode = 0,
                                      scale = "estimate",
                                      g_init = NULL,
                                      fix_g = FALSE,
                                      output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_ash <- function(x,
                     s = 1,
                     mode = 0,
                     scale = "estimate",
                     g_init = NULL,
                     fix_g = FALSE,
                     output = output_default(),
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

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_npmle <- function(x,
                       s = 1,
                       scale = "estimate",
                       g_init = NULL,
                       fix_g = FALSE,
                       output = output_default(),
                       optmethod = NULL,
                       control = NULL,
                       ...) {
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
                        call = match.call(),
                        ...))
}

#' @rdname ebnm_prior_families
#'
#' @export
#'
ebnm_deconvolver <- function(x,
                             s = 1,
                             scale = "estimate",
                             g_init = NULL,
                             fix_g = FALSE,
                             output = output_default(),
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
