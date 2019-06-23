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
                          ...) {
  return(ebnm_ash_workhorse(x = x,
                            s = s,
                            mode = mode,
                            scale = scale,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            call = match.call(),
                            mixcompdist = "halfuniform",
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
                                    ...) {
  return(ebnm_ash_workhorse(x = x,
                            s = s,
                            mode = mode,
                            scale = scale,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            call = match.call(),
                            mixcompdist = "uniform",
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
                                      ...) {
  return(ebnm_ash_workhorse(x = x,
                            s = s,
                            mode = mode,
                            scale = scale,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            call = match.call(),
                            mixcompdist = "+uniform",
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
                                      ...) {
  return(ebnm_ash_workhorse(x = x,
                            s = s,
                            mode = mode,
                            scale = scale,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            call = match.call(),
                            mixcompdist = "-uniform",
                            ...))
}
