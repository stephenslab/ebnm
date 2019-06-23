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
                                      ...) {
  return(ebnm_ash_workhorse(x = x,
                            s = s,
                            mode = mode,
                            scale = scale,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            call = match.call(),
                            mixcompdist = "normal",
                            ...))
}
