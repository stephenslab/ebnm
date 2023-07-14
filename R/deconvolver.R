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

  if (!is.null(call$mode)) {
    warning("mode parameter is ignored by ebnm_deconvolver.")
    call$mode <- NULL
  }

  if (is.null(g_init)) {
    g_init <- init_g_for_npmle(x, s, scale, force_pointmass = TRUE)
  }

  if (!fix_g) {
    tau_grid <- g_init$mean
    deconv_res <- do.call(deconvolveR::deconv,
                          c(list(tau = tau_grid,
                                 X = x,
                                 family = "Normal"),
                            control,
                            ...))
    pi <- deconv_res$stats[, "g"]
    g_init <- ashr::normalmix(pi = pi, mean = tau_grid, sd = 0)
  }

  orig_scale <- call$scale
  call$scale <- NULL
  ebnm_res <- ebnm_workhorse(x = x,
                             s = s,
                             mode = 0,
                             scale = "estimate",
                             g_init = g_init,
                             fix_g = TRUE,
                             output = output,
                             optmethod = NULL,
                             control = NULL,
                             prior_family = "npmle",
                             call = call)
  call$scale <- orig_scale

  if (!fix_g) {
    attr(ebnm_res$log_likelihood, "df") <- length(pi) - 1
  }

  return(ebnm_res)
}
