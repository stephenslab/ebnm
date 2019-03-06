#' @describeIn ebnm Solve the EBNM problem using a point-normal prior.
#'
#' @export
#'
ebnm_point_normal <- function (x,
                               s = 1,
                               g = list(mu = 0),
                               fixg = FALSE,
                               fix_pi0 = FALSE,
                               fix_mu = TRUE,
                               norm = NULL,
                               control = NULL,
                               output = NULL) {
  output <- set_output(output)

  if (fixg && is.null(g)) {
    stop("Must specify g if fixg = TRUE")
  }
  if (fix_pi0 && is.null(g$pi0)) {
    stop("Must specify g$pi0 if fix_pi0 = TRUE")
  }
  if (!is.null(g$a) && (g$a <= 0)) {
    stop("Invalid choice of g$a")
  }
  if (!is.null(g$pi0) && (g$pi0 < 0 || g$pi0 > 1)) {
    stop("Invalid choice of g$pi0")
  }

  if (all(is.infinite(s))) {
    stop("Impossible to fit g when all standard errors are infinite")
  }

  # Scale for stability, but need to be careful with log-likelihood:
  if (is.null(norm)) {
    pos_idx <- (is.finite(s) & s > 0)
    if (sum(pos_idx) == 0) {
      norm <- 1
    } else {
      norm <- mean(s[pos_idx])
    }
  }

  s <- s / norm
  x <- x / norm
  if (!is.null(g) && !is.null(g$mu)) {
    g$mu <- g$mu / norm
  }
  if (!is.null(g) && !is.null(g$a)) {
    g$a <- g$a * norm^2
  }

  # Don't use data with infinite SEs when estimating g:
  if (any(is.infinite(s))) {
    x_subset <- x[is.finite(s)]
    s_subset <- s[is.finite(s)]
  } else {
    x_subset <- x
    s_subset <- s
  }

  if (!fixg) {
    if (fix_pi0 & fix_mu) {
      g <- mle_point_normal_logscale_fixed_pi0_and_mu(x_subset, s_subset, g, control)
    } else if (fix_pi0 & !fix_mu) {
      g <- mle_point_normal_logscale_fixed_pi0(x_subset, s_subset, g, control)
    } else if (!fix_pi0 & fix_mu) {
      g <- mle_point_normal_logscale_fixed_mu(x_subset, s_subset, g, control)
    } else {
      g <- mle_point_normal_logscale_grad(x_subset, s_subset, g, control)
    }
  }

  w <- 1 - g$pi0
  a <- g$a
  mu <- g$mu

  # Compute return values, taking care to adjust back to original scale:
  retlist <- list()

  if ("result" %in% output) {
    result <- compute_summary_results_point_normal(x, s, w, a, mu)
    result$PosteriorMean <- result$PosteriorMean * norm
    result$PosteriorMean2 <- result$PosteriorMean2 * norm^2
    retlist <- c(retlist, list(result = result))
  }

  if ("fitted_g" %in% output) {
    fitted_g <- list(pi0 = g$pi0, a = g$a / norm^2, mu = g$mu * norm)
    retlist <- c(retlist, list(fitted_g = fitted_g))
  }

  if ("loglik" %in% output) {
    loglik <- loglik_point_normal(x_subset, s_subset, w, a, mu)
    loglik <- loglik - length(x_subset) * log(norm)
    retlist <- c(retlist, list(loglik = loglik))
  }

  if ("post_sampler" %in% output) {
    retlist <- c(retlist, list(post_sampler = function(nsamp) {
      post_sampler_point_normal(x, s, w, a, mu, nsamp) * norm
    }))
  }

  return(retlist)
}
