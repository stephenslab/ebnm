#' @importFrom ashr normalmix
#'
rebayes_workhorse <- function(x = x,
                              s = s,
                              mode = mode,
                              scale = scale,
                              g_init = g_init,
                              fix_g = fix_g,
                              output = output,
                              control = control,
                              call = call,
                              ...) {
  if (!requireNamespace("REBayes", quietly = TRUE)) {
    stop("Package REBayes must be installed to use optmethod = 'REBayes'.")
  }

  if (!is.null(call$mode)) {
    warning("mode parameter is ignored by ebnm_npmle.")
  }

  if (!fix_g && !is.null(g_init)) {
    stop("If optmethod = 'REBayes', then g_init must be NULL.")
  }

  if ("v" %in% names(list(...))) {
    stop("REBayes parameter 'v' is not supported. Please use 'scale' instead.")
  }

  if (fix_g) {
    df <- 0
  } else {
    if (scale == "estimate") {
      g <- init_g_for_npmle(x, s, force_pointmass = TRUE)
      n_gridpts <- length(g$pi)
    } else {
      n_gridpts <- (max(x) - min(x)) / scale
    }

    rebayes_res <- do.call(REBayes::GLmix,
                           c(list(x = x,
                                  sigma = s,
                                  v = n_gridpts),
                             control,
                             ...))

    g_init <- ashr::normalmix(pi = rebayes_res$y, mean = rebayes_res$x, sd = 0)
    df <- length(g_init$pi) - 1
    fix_g <- TRUE
    call$scale <- NULL
  }

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

  if (df > 0) {
    ebnm_res <- add_llik_to_retlist(ebnm_res, ebnm_res[[llik_ret_str()]], x, df = df)
  }

  return(ebnm_res)
}
