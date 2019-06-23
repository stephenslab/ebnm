#' @describeIn ebnm Solves the EBNM problem using a point-laplace prior.
#'
#' @export
#'
ebnm_point_laplace <- function (x,
                                s = 1,
                                mode = 0,
                                scale = "estimate",
                                g_init = NULL,
                                fix_g = FALSE,
                                output = output_default()) {

  if (mode != 0) {
    stop("Option to estimate mode not yet implemented for 'point_laplace' ",
         "priors.")
  }
  if (!identical(pmatch(scale, "estimate"), 1L)) {
    stop("Option to fix scale not yet implemented for 'point_laplace' priors.")
  }
  if (!is.null(g_init) && !fix_g) {
    stop("Option to intialize from g not yet implemented for 'point_laplace' ",
         "priors.")
  }
  if ("post_sampler" %in% output) {
    stop("Posterior sampler not yet implemented for 'point_laplace' priors.")
  }
  if (any(is.infinite(s))) {
    stop("Infinite SEs not yet implemented for 'point_laplace' priors.")
  }
  if (any(s == 0)) {
    stop("Zero SEs not yet implemented for 'point_laplace' priors.")
  }

  # TODO: could consider making more stable this way? But might have to be
  #   careful with log-likelihood.
  # m_sdev <- mean(s)
  # s <- s / m_sdev
  # x <- x / m_sdev

  # Estimate g from data
  if (!fix_g) {
    g <- mle_point_laplace(x, s)
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

	# Compute return values
	if ("result" %in% output) {
	  result <- summary_results_point_laplace(x, s, w, a)
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
	    loglik <- -g$val
	  }
	  retlist <- c(retlist, list(loglik = loglik))
	}

	return(retlist)
}

# Constructor for laplacemix class.
laplacemix <- function(pi, mean, scale) {
  structure(data.frame(pi, mean, scale), class="laplacemix")
}
