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
                                      call,
                                      ...) {
  if (length(setdiff(names(list(...)), "gridmult")) > 0) {
    warning("All additional parameters other than 'gridmult' are ignored for ",
            "scale mixtures of normal prior families. To use 'ashr' parameters, ",
            "use function 'ebnm_ash' with 'mixcompdist = \"normal\"'.")
  }

  if (!is.null(g_init)) {
    if (!inherits(g_init, "normalmix")) {
      stop("g_init must be NULL or an object of class normalmix.")
    }
    if (!is.null(call$mode) || !is.null(call$scale)) {
      warning("mode and scale parameters are ignored when g_init is supplied.")
    }
    min_s <- min(s)
    mode <- g_init$mean
    if (max(mode) - min(mode) < 1e-6 * min_s) {
      mode <- mode[1]
    }
    scale <- g_init$sd
    if (max(scale) - min(scale) < 1e-6 * min_s) {
      scale <- scale[1]
    }
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
                                            call = NULL,
                                            ...)
      return(ebnm_res[[llik_ret_str()]])
    }
    mode_opt_res <- optimize(mode_opt_fn, c(min(x), max(x)), maximum = TRUE)
    mode <- mode_opt_res$maximum
  }

  if (identical(scale, "estimate")) {
    if ("gridmult" %in% names(list(...))) {
      gridmult <- list(...)$gridmult
      scale <- get_ashr_grid(x, s, mode, gridmult)
    } else {
      scale <- ebnm_scale_normalmix(x, s, mode)
    }
  }

  n_mixcomp <- max(length(mode), length(scale))
  n_obs     <- length(x)

  # Adapted from ashr:::estimate_mixprop.
  if (length(s) == 1) {
    s <- rep(s, n_obs)
  }

  if (length(scale) > 1) {
    sigma2 <- outer(s^2, scale^2, `+`)
  } else {
    sigma2 <- s^2 + scale^2
  }
  if (length(mode) > 1) {
    llik_mat <- -0.5 * (log(sigma2) + outer(x, -mode, `+`)^2 / sigma2)
  } else {
    llik_mat <- -0.5 * (log(sigma2) + (x - mode)^2 / sigma2)
  }

  if (length(scale) == 1 && length(mode) == 1) {
    llik_mat <- matrix(llik_mat, ncol = 1)
  }

  llik_norms <- apply(llik_mat, 1, max)

  L_mat <- exp(llik_mat - llik_norms)

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
                          mean = rep(mode, length.out = n_mixcomp),
                          sd = rep(scale, length.out = n_mixcomp))
  }

  # Compute results.
  retlist <- list()

  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, x, s)
  }

  if (posterior_in_output(output) || sampler_in_output(output)) {
    posterior <- list()

    comp_postprob  <- L_mat * matrix(pi_est, n_obs, n_mixcomp, byrow = TRUE)
    comp_postprob  <- comp_postprob / rowSums(comp_postprob)

    if (length(mode) > 1) {
      comp_postmean <- outer(s^2, mode)
    } else {
      comp_postmean <- mode * s^2
    }
    if (length(scale) > 1) {
      comp_postmean  <- (comp_postmean + outer(x, scale^2)) / sigma2
      comp_postmean2 <- comp_postmean^2 + outer(s^2, scale^2) / sigma2
    } else {
      comp_postmean  <- (comp_postmean + x * scale^2) / sigma2
      comp_postmean2 <- comp_postmean^2 + s^2 * scale^2 / sigma2
    }

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
      } else {
        probzero <- 0
      }
      probneg <- 1 - (probpos + probzero)
      posterior$lfsr <- probzero + pmin(probneg, probpos)
    }

    retlist <- add_posterior_to_retlist(retlist, posterior, output, x)
  }

  if (g_in_output(output)) {
    retlist <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    loglik  <- sum(log(L_mat %*% pi_est))
    loglik  <- loglik + sum(llik_norms) - n_obs * log(2 * pi) / 2
    df      <- (1 - fix_g) * (length(fitted_g$pi) - 1)
    retlist <- add_llik_to_retlist(retlist, loglik, x, df = df)
  }

  if (sampler_in_output(output)) {
    # Adapted from ashr::post_sample.normalmix.
    post_sampler <- function(nsamp) {
      mixcomp <- apply(comp_postprob, 1, function(prob) {
        sample(1:n_mixcomp, nsamp, replace = TRUE, prob = prob)
      })
      idx <- rep(1:n_obs, each = nsamp) + (mixcomp - 1) * (n_obs)
      samp_means <- comp_postmean[idx]
      samp_sds <- sqrt(comp_postmean2[idx] - samp_means^2)
      samp <- rnorm(nsamp * n_obs, mean = samp_means, sd = samp_sds)
      samp <- matrix(samp, nrow = nsamp, ncol = n_obs)
      colnames(samp) <- names(x)
      return(samp)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(retlist)
}
