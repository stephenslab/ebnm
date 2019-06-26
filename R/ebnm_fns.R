#' @describeIn ebnm Solves the EBNM problem using a point-normal prior.
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
                              control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        control = control,
                        prior_family = "point_normal",
                        call = match.call()))
}

#' @describeIn ebnm Solves the EBNM problem using a point-Laplace prior.
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
                               control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        control = control,
                        prior_family = "point_laplace",
                        call = match.call()))
}

#' @describeIn ebnm Solves the EBNM problem using a normal prior (with no
#'   point mass).
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
                        control = NULL) {
  return(ebnm_workhorse(x = x,
                        s = s,
                        mode = mode,
                        scale = scale,
                        g_init = g_init,
                        fix_g = fix_g,
                        output = output,
                        control = control,
                        prior_family = "normal",
                        call = match.call()))
}

#' @describeIn ebnm Solves the EBNM problem using a scale mixture of normals.
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
                        control = control,
                        prior_family = "normal_scale_mixture",
                        call = match.call()))
}

#' @describeIn ebnm Solves the EBNM problem using a unimodal distribution.
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
                        control = control,
                        prior_family = "unimodal",
                        call = match.call(),
                        ...))
}

#' @describeIn ebnm Solves the EBNM problem using a symmetric unimodal
#'   distribution.
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
                        control = control,
                        prior_family = "unimodal_symmetric",
                        call = match.call(),
                        ...))
}

#' @describeIn ebnm Solves the EBNM problem using a unimodal distribution with
#'   support constrained to be greater than the mode.
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
                        control = control,
                        prior_family = "unimodal_nonnegative",
                        call = match.call(),
                        ...))
}

#' @describeIn ebnm Solves the EBNM problem using a unimodal distribution with
#'   support constrained to be less than the mode.
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
                        control = control,
                        prior_family = "unimodal_nonpositive",
                        call = match.call(),
                        ...))
}

#' @describeIn ebnm A generic function for solving the EBNM problem using an
#'   adaptive shrinkage (ash) prior.
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
                        control = control,
                        prior_family = "ash",
                        call = match.call(),
                        ...))
}
