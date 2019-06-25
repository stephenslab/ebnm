

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
                              pointmass) {
  if (!is.null(g_init)) {
    if (!inherits(g_init, "normalmix")) {
      stop("g_init must be NULL or an object of class ashr::normalmix.")
    }
    ncomp <- length(g_init$pi)
    if ((ncomp == 0) || (ncomp > 2)) {
      stop("g_init does not have the correct number of components.")
    }
    g <- list(pi0 = g_init$pi[1],
              a = 1 / g_init$sd[ncomp]^2,
              mu = g_init$mean[1])
  } else {
    g <- list()
  }

  if (pointmass) {
    fix_pi0 <- fix_g
  } else {
    fix_pi0 <- TRUE
    g$pi0 <- 0
  }

  # Allow partial matching for mode and scale.
  if (identical(pmatch(mode, "estimate"), 1L)) {
    fix_mu <- fix_g
  } else if (is.numeric(mode) && (length(mode) == 1)) {
    # Use all.equal to allow for numerical error.
    if (!is.null(g$mu) && !isTRUE(all.equal(g$mu, mode))) {
      stop("If mode and g_init$mean are both supplied, they must agree.")
    }
    fix_mu <- TRUE
    g$mu <- mode
  } else {
    stop("Invalid argument to mode.")
  }

  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }

  if (identical(pmatch(scale, "estimate"), 1L)) {
    fix_a <- fix_g
  } else if (is.numeric(scale) && (length(scale) == 1) && (scale > 0)) {
    if (!is.null(g$a) && !isTRUE(all.equal(g$a, 1 / scale^2))) {
      stop("If scale and g_init$sd are both supplied, they must agree.")
    }
    fix_a <- TRUE
    g$a <- 1 / scale^2
  } else {
    stop("Invalid argument to scale.")
  }

  if (fix_pi0 && fix_mu && fix_a) {
    fix_g <- TRUE
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
      g <- mle_point_normal(x_optset, s_optset, g, control,
                            fix_pi0, fix_a, fix_mu)
    }
  }

  pi0 <- g$pi0
  w <- 1 - pi0
  a <- g$a
  mu <- g$mu

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
