#' @describeIn ebnm Solves the EBNM problem using a scale mixture of normals.
#'
#' @importFrom ashr normalmix
#'
#' @export
#'
ebnm_normal_scale_mixture <- function(x,
                                      s = 1,
                                      mode = 0,
                                      scale = "estimate",
                                      g_init = NULL,
                                      fix_g = FALSE,
                                      output = output_default(),
                                      control = list(),
                                      pointmass = TRUE,
                                      grid_mult = sqrt(2),
                                      ...) {
  # If mode is "estimate" then just call into ashr.
  if (is.null(g_init) && identical(pmatch(mode, "estimate"), 1L)) {
    return(ebnm_ash_workhorse(x = x,
                              s = s,
                              mode = mode,
                              scale = scale,
                              fixg = fixg,
                              output = output,
                              call = match.call(),
                              mixcompdist = "normal",
                              pointmass = pointmass,
                              mult = grid_mult,
                              ...))
  }

  check_args(x, s, g_init, fix_g, output)

  if (!is.null(g_init)) {
    if (!inherits(g_init, "normalmix")) {
      stop("g_init must be NULL or an object of class ashr::normalmix.")
    }
    if (!identical(pmatch(scale, "estimate"), 1L)) {
      warning("An initial g was supplied. Ignoring scale parameter.")
    }
    scale <- g_init$sd
  } else {
    # Adapted from ashr:::autoselect.mixsd.
    sigmamin <- min(s[s > 0]) / 10
    sigmamax <- max(8 * sigmamin, 2 * sqrt(max(x^2 - s^2, 0)))
    npoint <- ceiling(log2(sigmamax / sigmamin) / log2(grid_mult))
    scale <- grid_mult^((-npoint):0) * sigmamax
    if (pointmass) {
      scale <- c(0, scale)
    }
  }

  n_mixcomp <- length(scale)
  n_obs     <- length(x)

  # Adapted from ashr:::estimate_mixprop. Rows correspond to observations;
  #   columns, to mixture components.
  if (length(s) == 1) {
    sigmavec <- s^2 + scale^2
    llik_mat <- -0.5 * t(log(sigmavec) + outer(1 / sigmavec, (x - mode)^2))
  } else {
    sigmamat <- outer(s^2, scale^2, `+`)
    llik_mat <- -0.5 * (log(sigmamat) + (x - mode)^2 / sigmamat)
  }
  llik_norms <- apply(llik_mat, 1, max)
  llik_mat   <- llik_mat - llik_norms

  if (fix_g) {
    fitted_g <- g_init
    est_pi <- g_init$pi
  } else {
    if (is.null(g_init)) {
      x0 <- rep(1, n_mixcomp)
    } else {
      x0 <- g_init$pi
    }

    L <- exp(llik_mat)
    nonzero_cols <- (apply(L, 2, max) > 0)
    if (!all(nonzero_cols)) {
      x0 <- x0[nonzero_cols]
      L  <- L[, nonzero_cols, drop = FALSE]
    }

    control0 <- list(verbose = FALSE)
    control  <- modifyList(control0, control, keep.null = TRUE)
    optres   <- mixsqp::mixsqp(L, x0 = x0, control = control)

    est_pi <- rep(0, n_mixcomp)
    est_pi[nonzero_cols] <- pmax(optres$x, 0)
    est_pi <- est_pi / sum(est_pi)
    fitted_g <- normalmix(pi = est_pi,
                          mean = rep(mode, n_mixcomp),
                          sd = scale)
  }

  retlist <- list()

  if (any(c("result", "lfsr", "post_sampler") %in% output)) {
    comp_postprob <- t(exp(t(llik_mat) + log(est_pi)))
    comp_postprob <- comp_postprob / rowSums(comp_postprob)
    if (length(s) == 1) {
      comp_postmean  <- mode * s^2 + t(outer(scale^2, x) / sigmavec)
      comp_postmean2 <- t(t(comp_postmean)^2 + s^2 * scale^2 / sigmavec)
    } else {
      comp_postmean  <- mode * s^2 + outer(x, scale^2) / sigmamat
      comp_postmean2 <- comp_postmean^2 + outer(s^2, scale^2) / sigmamat
    }

    result <- list()
    if ("result" %in% output) {
      result$posterior_mean  <- rowSums(comp_postprob * comp_postmean)
      result$posterior_mean2 <- rowSums(comp_postprob * comp_postmean2)
    }
    if ("lfsr" %in% output) {
      comp_probpos <- pnorm(0,
                            mean = comp_postmean,
                            sd = sqrt(comp_postmean2 - comp_postmean^2),
                            lower.tail = FALSE)
      probpos <- rowSums(comp_postprob * comp_probpos)
      zero_comps <- which(fitted_g$sd == 0)
      if (length(zero_comps) > 0) {
        probzero <- rowSums(comp_postprob[, zero_comps, drop = FALSE])
      }
      probneg <- 1 - (probpos + probzero)
      result$lfsr <- probzero + pmin(probneg, probpos)
    }
    retlist <- c(retlist, list(result = data.frame(result)))
  }

  if ("fitted_g" %in% output) {
    retlist <- c(retlist, list(fitted_g = fitted_g))
  }

  if ("loglik" %in% output) {
    loglik <- sum(log(exp(llik_mat) %*% est_pi))
    loglik <- loglik + sum(llik_norms) - n_obs * log(2 * pi) / 2
    retlist <- c(retlist, list(loglik = loglik))
  }

  if ("post_sampler" %in% output) {
    # Adapted from ashr::post_sample.normalmix.
    post_sampler <- function(nsamp) {
      mixcomp <- apply(comp_postprob, 1, function(prob) {
        sample(1:n_mixcomp, nsamp, replace = TRUE, prob = prob)
      })
      idx <- rep(1:n_obs, each = nsamp) + (mixcomp - 1) * (n_obs)
      samp_means <- comp_postmean[idx]
      samp_sds   <- sqrt(comp_postmean2[idx] - samp_means^2)
      samp       <- rnorm(nsamp * n_obs, mean = samp_means, sd = samp_sds)
      return(matrix(samp, nrow = nsamp, ncol = n_obs))
    }
    retlist <- c(retlist, list(post_sampler = post_sampler))
  }

  return(retlist)
}
