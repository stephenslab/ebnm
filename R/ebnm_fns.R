#' @rdname ebnm_prior_families
#' 
#' @title Solve the EBNM Problem Using a Specified Prior
#'
#' @details \code{ebnm_point_normal} uses a point-normal prior;
#'   \code{ebnm_point_laplace} uses a point-Laplace prior;
#'   \code{ebnm_point_exponential} uses a point-exponential prior;
#'   \code{ebnm_normal} uses a normal prior; \code{ebnm_horseshoe} uses
#'   a horseshoe prior; \code{ebnm_normal_scale_mixture} uses a scale
#'   mixture of normals; \code{ebnm_unimodal} uses a unimodal
#'   distribution; \code{ebnm_unimodal_symmetric} uses a symmetric
#'   unimodal distribution; \code{ebnm_unimodal_nonnegative} uses a
#'   nonnegative unimodal distribution; \code{ebnm_unimodal_nonpositive}
#'   uses a nonpositive unimodal distribution; \code{ebnm_ash} uses an
#'   adaptive shrinkage (\dQuote{ash}) prior; \code{ebnm_npmle} solves
#'   the EBNM problem using a completely nonparametric approach; and
#'   \code{ebnm_deconvoler} solves the EBNM problem using the
#'   \code{deconvolveR} pacakage.
#'
#' @inheritParams ebnm
#'
#' @return Describe return value here.
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
