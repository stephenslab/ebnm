ebnm_pl_workhorse <- function(x,
                              s,
                              mode,
                              scale,
                              g_init,
                              fix_g,
                              output,
                              optmethod,
                              use_grad,
                              use_hess,
                              control,
                              call) {
  if (length(scale) != 1) {
    stop("Argument 'scale' must be either 'estimate' or a scalar.")
  }

  check_g_init(g_init,
               fix_g,
               mode = mode,
               scale = scale,
               pointmass = TRUE,
               call = call,
               class_name = "laplacemix",
               scale_name = "scale")

  fix_pi0 <- FALSE
  fix_a   <- !identical(scale, "estimate")
  fix_mu  <- !identical(mode, "estimate")

  if (!is.null(g_init) && length(g_init$pi) == 1) {
    g <- list(pi0 = 0,
              a = 1 / g_init$scale,
              mu = g_init$mean)
  } else if (!is.null(g_init) && length(g_init$pi) == 2) {
    g <- list(pi0 = g_init$pi[1],
              a = 1 / g_init$scale[2],
              mu = g_init$mean[1])
  } else {
    g <- list()
    if (fix_pi0) {
      g$pi0 <- 0
    }
    if (fix_a) {
      g$a <- 1 / scale
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
    if (fix_a && fix_mu) {
      g <- mle_point_laplace_fixa(x_optset, s_optset, g, control)
    } else {
      g <- mle_point_laplace(x_optset, s_optset, g, control,
                             fix_pi0, fix_a, fix_mu,
                             optmethod, use_grad, use_hess)
    }
  }

  pi0 <- g$pi0
  w   <- 1 - g$pi0
  a   <- g$a
  mu  <- g$mu

  retlist <- list()

  if (posterior_in_output(output)) {
    posterior <- summary_results_point_laplace(x, s, w, a, mu, output)
    retlist   <- add_posterior_to_retlist(retlist, posterior, output)
  }

  if (g_in_output(output)) {
    if (pi0 == 0) {
      fitted_g <- laplacemix(pi = 1, mean = mu, scale = 1 / a)
    } else {
      fitted_g <- laplacemix(pi = c(pi0, w),
                             mean = rep(mu, 2),
                             scale = c(0, 1 / a))
    }
    retlist <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    if (fix_g) {
      loglik <- loglik_point_laplace(x_optset, s_optset, w, a, mu)
    } else {
      loglik <- g$val
    }
    retlist <- add_llik_to_retlist(retlist, loglik)
  }

  if (sampler_in_output(output)) {
    post_sampler <- function(nsamp) {
      post_sampler_point_laplace(x, s, w, a, mu, nsamp)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(retlist)
}

#' Constructor for laplacemix class
#'
#' Creates a finite mixture of Laplace distributions.
#'
#' @param pi A vector of mixture proportions.
#' @param mean A vector of means.
#' @param scale A vector of scale parameters.
#'
#' @export
#'
laplacemix <- function(pi, mean, scale) {
  structure(data.frame(pi, mean, scale), class="laplacemix")
}
