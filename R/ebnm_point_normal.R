#' @describeIn ebnm Solve the EBNM problem using a point-normal prior.
#'
#' @export
#'
ebnm_point_normal <- function (x,
                               s = 1,
                               g = list(),
                               fixg = FALSE,
                               fix_pi0 = FALSE,
                               fix_mu = TRUE,
                               norm = NULL,
                               control = NULL,
                               output = NULL,
                               use_cpp = FALSE) {
  output <- set_output(output)
  check_args(x, s, g, fixg, output)

  # If mu is fixed but unspecified, fix it at zero.
  if ((fixg || fix_mu) && is.null(g$mu)) {
    g$mu <- 0
  }
  if (fix_pi0 && is.null(g$pi0)) {
    stop("Must specify g$pi0 if fix_pi0 = TRUE.")
  }

  x_optset <- x
  s_optset <- s
  # Don't use data with infinite SEs when estimating g.
  if (any(is.infinite(s))) {
    x_optset <- x[is.finite(s)]
    s_optset <- s[is.finite(s)]
  }

  if (!fixg) {
    if (use_cpp) {
      g <- cpp_mle_point_normal(x_optset, s_optset, g, control, fix_pi0, fix_mu)
    } else {
      g <- mle_point_normal_logscale_grad(x_optset, s_optset, g, control,
                                          fix_pi0, fix_mu)
    }
  }

  w <- 1 - g$pi0
  a <- g$a
  mu <- g$mu

  retlist <- list()
  if ("result" %in% output) {
    result <- compute_summary_results_point_normal(x, s, w, a, mu)
    retlist <- c(retlist, list(result = result))
  }
  if ("fitted_g" %in% output) {
    fitted_g <- list(pi0 = g$pi0, a = g$a, mu = g$mu)
    retlist <- c(retlist, list(fitted_g = fitted_g))
  }
  if ("loglik" %in% output) {
    if (!fixg) {
      loglik <- -g$val
    } else {
      loglik <- loglik_point_normal(x_optset, s_optset, w, a, mu)
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
