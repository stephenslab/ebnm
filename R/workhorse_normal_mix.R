# This does basically the same thing as ashr::ash, but is re-organized to
#   avoid doing any given matrix computation more than a single time.
#
#' @importFrom ashr normalmix
#' @importFrom mixsqp mixsqp
#' @importFrom stats optimize rnorm pnorm
#' @importFrom utils modifyList
#'
ebnm_normal_mix_workhorse <- function(x,
                                      s,
                                      mode,
                                      scale,
                                      g_init,
                                      fix_g,
                                      output,
                                      control,
                                      pointmass,
                                      grid_mult,
                                      call) {
  if (!is.null(g_init)) {
    if (!inherits(g_init, "normalmix")) {
      stop("g_init must be NULL or an object of class normalmix.")
    }
    if (!is.null(call$mode) || !is.null(call$scale)) {
      warning("mode and scale parameters are ignored when g_init is supplied.")
    }
    mode  <- g_init$mean[1]
    scale <- g_init$sd
  }

  if (identical(mode, "estimate")) {
    # Adapted from ash.estmode.
    mode_opt_fn <- function(m) {
      ebnm_res <- ebnm_normal_mix_workhorse(x = x,
                                            s = s,
                                            mode = m,
                                            scale = scale,
                                            g_init = NULL,
                                            fix_g = FALSE,
                                            output = llik_arg_str(),
                                            control = control,
                                            pointmass = pointmass,
                                            grid_mult = grid_mult,
                                            call = NULL)
      return(ebnm_res[[llik_ret_str()]])
    }
    mode_opt_res <- optimize(mode_opt_fn, c(min(x), max(x)), maximum = TRUE)
    mode <- mode_opt_res$maximum
  }

  if (identical(scale, "estimate")) {
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

  # Adapted from ashr:::estimate_mixprop.
  if (length(s) == 1) {
    s <- rep(s, n_obs)
  }
  sigmamat   <- outer(s^2, scale^2, `+`)
  llik_mat   <- -0.5 * (log(sigmamat) + (x - mode)^2 / sigmamat)
  llik_norms <- apply(llik_mat, 1, max)
  L_mat      <- exp(llik_mat - llik_norms)

  if (fix_g) {
    fitted_g <- g_init
    pi_est   <- g_init$pi
  } else {
    if (is.null(g_init)) {
      pi_init <- rep(1, n_mixcomp)
    } else {
      pi_init <- g_init$pi
    }

    nonzero_cols <- (apply(L_mat, 2, max) > 0)
    if (all(nonzero_cols)) {
      x0 <- pi_init
      L  <- L_mat
    } else {
      x0 <- pi_init[nonzero_cols]
      L  <- L_mat[, nonzero_cols, drop = FALSE]
    }

    control0 <- list(verbose = FALSE)
    control  <- modifyList(control0, control, keep.null = TRUE)
    optres   <- mixsqp(L = L, x0 = x0, control = control)

    pi_est <- rep(0, n_mixcomp)
    pi_est[nonzero_cols] <- pmax(optres$x, 0)
    pi_est <- pi_est / sum(pi_est)
    fitted_g <- normalmix(pi = pi_est,
                          mean = rep(mode, n_mixcomp),
                          sd = scale)
  }

  # Compute results.
  retlist <- list()

  if (posterior_in_output(output)) {
    posterior <- list()

    comp_postprob  <- L_mat * matrix(pi_est, n_obs, n_mixcomp, byrow = TRUE)
    comp_postprob  <- comp_postprob / rowSums(comp_postprob)
    comp_postmean  <- mode * s^2 + outer(x, scale^2) / sigmamat
    comp_postmean2 <- comp_postmean^2 + outer(s^2, scale^2) / sigmamat

    if (result_in_output(output)) {
      posterior$mean  <- rowSums(comp_postprob * comp_postmean)
      posterior$mean2 <- rowSums(comp_postprob * comp_postmean2)
      posterior$sd    <- sqrt(pmax(0, posterior$mean2 - posterior$mean^2))
    }

    if (lfsr_in_output(output)) {
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
      posterior$lfsr <- probzero + pmin(probneg, probpos)
    }

    retlist <- add_posterior_to_retlist(retlist, posterior, output)
  }

  if (g_in_output(output)) {
    retlist <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    loglik <- sum(log(L_mat %*% pi_est))
    loglik <- loglik + sum(llik_norms) - n_obs * log(2 * pi) / 2
    retlist <- add_llik_to_retlist(retlist, loglik)
  }

  if (sampler_in_output(output)) {
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
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(retlist)
}
