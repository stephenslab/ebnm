#' @describeIn ebnm Solves the EBNM problem using a mixture of half-uniforms.
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

#' @describeIn ebnm Solves the EBNM problem using a mixture of uniforms.
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
                            ...))
}

#' @describeIn ebnm Solves the EBNM problem using a +uniform mixture.
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

#' @describeIn ebnm Solves the EBNM problem using a -uniform mixture.
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
