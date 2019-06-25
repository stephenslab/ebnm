ebnm_pl_workhorse <- function(x,
                              s,
                              mode,
                              scale,
                              g_init,
                              fix_g,
                              output,
                              control) {
  if (mode != 0) {
    stop("Option to estimate mode not yet implemented for 'point_laplace' ",
         "priors.")
  }

  if (any(s == 0)) {
    stop("Handling of SEs equal to zero not yet implemented for ",
         "'point_laplace' priors.")
  }

  if (!is.null(g_init)) {
    if (!inherits(g_init, "laplacemix")) {
      stop("g_init must be NULL or an object of class laplacemix.")
    }
    ncomp <- length(g_init$pi)
    if (ncomp != 2) {
      stop("g_init does not have the correct number of components.")
    }
    g <- list(pi0 = g_init$pi[1],
              a = 1 / g_init$scale[2])
  } else {
    g <- list()
  }

  # Allow partial matching for scale.
  if (identical(pmatch(scale, "estimate"), 1L)) {
    fix_a <- fix_g
  } else if (is.numeric(scale) && (length(scale) == 1) && (scale > 0)) {
    if (!is.null(g$a) && !isTRUE(all.equal(g$a, 1 / scale))) {
      stop("If scale and g_init$scale are both supplied, they must agree.")
    }
    fix_a <- TRUE
    g$a <- 1 / scale
  } else {
    stop("Invalid argument to scale.")
  }

  x_optset <- x
  s_optset <- s
  # Don't use observations with infinite SEs when estimating g.
  if (any(is.infinite(s))) {
    x_optset <- x[is.finite(s)]
    s_optset <- s[is.finite(s)]
  }

  # Estimate g.
  if (!fix_g) {
    if (fix_a) {
      g <- mle_point_laplace_fixa(x_optset, s_optset, g)
    } else {
      g <- mle_point_laplace(x_optset, s_optset, g, control)
    }
  } else {
    if (!inherits(g_init, "laplacemix")) {
      stop("g_init must be NULL or an object of class laplacemix.")
    }
    g <- list(pi0 = g_init$pi[1], a = 1 / g_init$scale[2])
  }

  pi0 <- g$pi0
  w <- 1 - g$pi0
  a <- g$a

  retlist <- list()

  if ("result" %in% output || "lfsr" %in% output) {
    result <- summary_results_point_laplace(x, s, w, a, output)
    retlist <- c(retlist, list(result = result))
  }

  if ("fitted_g" %in% output) {
    retlist <- c(retlist,
                 list(fitted_g = laplacemix(pi = c(pi0, w),
                                            mean = rep(0, 2),
                                            scale = c(0, 1 / a))))
  }

  if ("loglik" %in% output) {
    if (fix_g) {
      loglik <- loglik_point_laplace(x, s, w, a)
    } else {
      loglik <- g$val
    }
    retlist <- c(retlist, list(loglik = loglik))
  }

  if ("post_sampler" %in% output) {
    retlist <- c(retlist, list(post_sampler = function(nsamp) {
      post_sampler_point_laplace(x, s, w, a, nsamp)
    }))
  }

  return(retlist)
}

# Constructor for laplacemix class.
laplacemix <- function(pi, mean, scale) {
  structure(data.frame(pi, mean, scale), class="laplacemix")
}
