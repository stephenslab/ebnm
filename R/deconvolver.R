#' @importFrom deconvolveR deconv
#' @importFrom ashr normalmix
#'
deconvolver_workhorse <- function(x = x,
                                  s = s,
                                  mode = mode,
                                  scale = scale,
                                  g_init = g_init,
                                  fix_g = fix_g,
                                  output = output,
                                  control = control,
                                  call = call,
                                  ...) {
  # deconvolveR takes z-scores:
  if (!isTRUE(all.equal(s, 1))) {
    stop("deconvolveR takes z-scores rather than observations and standard",
         " errors. Please convert to z-scores and set parameter 's = 1'.")
  }

  if (is.null(g_init)) {
    if (!is.null(call$mode)) {
      warning("mode parameter is ignored by ebnm_npmle.")
      call$mode <- NULL
    }

    g_init <- init_g_for_npmle(x, s, scale, force_pointmass = TRUE)
    tau_grid <- g_init$mean
    deconv_res <- do.call(deconvolveR::deconv,
                          c(list(tau = tau_grid,
                                 X = x,
                                 family = "Normal"),
                            control,
                            ...))
    call$scale <- NULL

    pi <- deconv_res$stats[, "g"]
    g_init <- ashr::normalmix(pi = pi, mean = tau_grid, sd = 0)
    fix_g <- TRUE
  }

  ebnm_res <- ebnm_workhorse(x = x,
                             s = s,
                             mode = 0,
                             scale = "estimate",
                             g_init = g_init,
                             fix_g = fix_g,
                             output = output,
                             optmethod = NULL,
                             control = NULL,
                             prior_family = "npmle",
                             call = call)

  return(ebnm_res)
}
