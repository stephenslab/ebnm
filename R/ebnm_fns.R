#' Solve the EBNM problem using a point-normal prior
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a point-Laplace prior
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a point-exponential prior
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a normal prior
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a scale mixture of normals
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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
                                      control = NULL) {
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
                        call = match.call()))
}

#' Solve the EBNM problem using a unimodal distribution
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a symmetric unimodal distribution
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a nonnegative unimodal distribution
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a nonpositive unimodal distribution
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using an ash opior
#'
#' A generic function for solving the EBNM problem using an
#'   adaptive shrinkage (ash) prior. See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
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

#' Solve the EBNM problem using a completely nonparametric problem
#'
#' See \code{\link{ebnm}} for details.
#'
#' @inheritParams ebnm
#'
#' @export
#'
ebnm_npmle <- function(x,
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
                        control = control,
                        prior_family = "npmle",
                        call = match.call(),
                        ...))
}
