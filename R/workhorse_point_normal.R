# The workhorse function is used by both ebnm_point_normal and ebnm_normal.
#
#' @importFrom ashr normalmix
#'
ebnm_pn_workhorse <- function(x,
                              s,
                              mode,
                              scale,
                              g_init,
                              fix_g,
                              output,
                              control,
                              pointmass,
                              call) {
  if (length(scale) != 1) {
    stop("Argument 'scale' must be either 'estimate' or a scalar.")
  }

  check_g_init(g_init,
               fix_g,
               pointmass = pointmass,
               call = call,
               class_name = "normalmix",
               scale_name = "sd")

  fix_pi0 <- !pointmass
  fix_a   <- !identical(scale, "estimate")
  fix_mu  <- !identical(mode, "estimate")

  if (fix_pi0 && fix_mu && fix_a) {
    fix_g <- TRUE
  }

  if (!is.null(g_init) && length(g_init$pi) == 1) {
    g <- list(pi0 = 0,
              a = 1 / g_init$sd^2,
              mu = g_init$mean)
  } else if (!is.null(g_init) && length(g_init$pi) == 2) {
    g <- list(pi0 = g_init$pi[1],
              a = 1 / g_init$sd[2]^2,
              mu = g_init$mean[1])
  } else {
    g <- list()
    if (fix_pi0) {
      g$pi0 <- 0
    }
    if (fix_a) {
      g$a <- 1 / scale^2
    }
    if (fix_mu) {
      g$mu <- mode
    }
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
    if (fix_pi0 && g$pi0 == 1) {
      g <- mle_point_only(x_optset, s_optset, g, fix_a, fix_mu)
    } else if (fix_pi0 && g$pi0 == 0) {
      g <- mle_normal(x_optset, s_optset, g, fix_a, fix_mu)
    } else {
      g <- mle_point_normal(x_optset, s_optset, g, control, fix_pi0, fix_a, fix_mu)
    }
  }

  pi0 <- g$pi0
  w   <- 1 - pi0
  a   <- g$a
  mu  <- g$mu

  retlist <- list()

  if ("result" %in% output || "lfsr" %in% output) {
    result <- summary_results_point_normal(x, s, w, a, mu, output)
    retlist <- c(retlist, list(result = result))
  }

  if ("fitted_g" %in% output) {
    if (pi0 == 0) {
      fitted_g <- normalmix(pi = 1, mean = mu, sd = sqrt(1 / a))
    } else {
      fitted_g <- normalmix(pi = c(pi0, w),
                            mean = rep(mu, 2),
                            sd = c(0, sqrt(1 / a)))
    }
    retlist <- c(retlist, list(fitted_g = fitted_g))
  }

  if ("loglik" %in% output) {
    if (fix_g) {
      loglik <- loglik_point_normal(x_optset, s_optset, w, a, mu)
    } else {
      loglik <- g$val
    }
    retlist <- c(retlist, list(loglik = loglik))
  }

  if ("post_sampler" %in% output) {
    retlist <- c(retlist, list(post_sampler = function(nsamp) {
      post_sampler_point_normal(x, s, w, a, mu, nsamp)
    }))
  }

  return(retlist)
}

# Used by both ebnm_pn_workhorse and ebnm_pl_workhorse.
check_g_init <- function(g_init, fix_g, pointmass, call, class_name, scale_name) {
  if (!is.null(g_init)) {
    if (!inherits(g_init, class_name)) {
      stop("g_init must be NULL or an object of class ", class_name, ".")
    }

    ncomp <- length(g_init$pi)
    if (!(ncomp == 1 || (pointmass && ncomp == 2))) {
      stop("g_init does not have the correct number of components.")
    }
    if (ncomp == 2 && g_init[[scale_name]][1] != 0) {
      stop("The first component of g_init must be a point mass.")
    }

    if (fix_g && (!is.null(call$mode) || !is.null(call$scale))) {
      warning("mode and scale parameters are ignored when g is fixed.")
    }

    if (!fix_g) {
      # all.equal allows for numerical error:
      if (!is.null(call$mode)
          && !identical(call$mode, "estimate")
          && !isTRUE(all.equal(g_init$mean[1], call$mode))) {
        stop("If mode is fixed and g_init is supplied, they must agree.")
      }
      g_scale <- g_init[[scale_name]][ncomp]
      if (!is.null(call$scale)
          && !identical(call$scale, "estimate")
          && !isTRUE(all.equal(g_scale, scale))) {
        stop("If scale is fixed and g_init is supplied, they must agree.")
      }
    }
  }
}
