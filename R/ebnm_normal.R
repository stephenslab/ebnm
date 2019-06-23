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
                        output = output_default()) {
  return(ebnm_pn_workhorse(x = x,
                           s = s,
                           mode = mode,
                           scale = scale,
                           g_init = g_init,
                           fix_g = fix_g,
                           output = output,
                           control = NULL,
                           pointmass = FALSE))
}
