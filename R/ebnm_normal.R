#' @describeIn ebnm Solve the EBNM problem using a normal prior.
#'
#' @export
#'
ebnm_normal <- function(x,
                        s = 1,
                        g_init = NULL,
                        fix_g = FALSE,
                        output = output_default(),
                        mode = 0,
                        sd = "estimate") {
  return(ebnm_pn_workhorse(x = x,
                           s = s,
                           g_init = g_init,
                           fix_g = fix_g,
                           output = output,
                           mode = mode,
                           sd = sd,
                           control = NULL,
                           pointmass = FALSE))
}
