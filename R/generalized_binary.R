### Author: Yusha Liu (with edits by J. Willwerscheid).
###   mode corresponds to mode of truncated normal component.
###   scale corresponds to sigma / mu.

ebnm_generalized_binary_defaults <- function(x, s) {
  min_mu <- min(s) / 10
  if (all(x <= 0)) {
    mu_init <- min_mu
    max_mu <- max(s)
  } else {
    mu_init <- mean(x[x > 0])
    max_mu <- max(x - s)
  }

  return(list(
    maxiter = 50,
    tol = 1e-3,
    wlist = c(1e-5, 1),
    mu_init = mu_init,
    mu_range = c(min_mu, max_mu)
  ))
}

#' @importFrom ashr tnormalmix
#' @importFrom stats pnorm optim
#'
gb_workhorse <- function(x,
                         s,
                         mode = "estimate",
                         scale = 0.1,
                         g_init,
                         fix_g,
                         output,
                         control,
                         call,
                         ...) {
  args <- ebnm_generalized_binary_defaults(x, s)
  for (arg in names(list(...))) {
    if (!(arg %in% names(args))) {
      warning("Argument ", arg, " is not recognized by ebnm_generalized_binary.")
    } else {
      args[[arg]] <- list(...)[[arg]]
    }
  }

  maxiter <- args$maxiter
  tol <- args$tol
  wlist <- args$wlist
  mu_init <- args$mu_init
  mu_range <- args$mu_range

  if (!is.null(g_init) && !inherits(g_init, "tnormalmix")) {
    stop("g_init must be NULL or an object of class tnormalmix")
  }

  if (!is.null(g_init)) {
    if (!(length(g_init$pi) == 2 &&
          g_init$mean[1] == 0 &&
          g_init$sd[1] == 0 &&
          g_init$a[2] == 0 &&
          g_init$b[2] == Inf)) {
      stop("g_init is a tnormalmix object, but it is not a generalized ",
           "binary prior (i.e., a mixture of a point mass at zero and a ",
           "truncated normal component with lower bound zero).")
    }
    if (!is.null(call$mode) || !is.null(call$scale)) {
      warning("mode and scale parameters are ignored when g_init is supplied.")
    }
    scale <- g_init$sd[2] / g_init$mean[2]
  } else if (identical(scale, "estimate")) {
    stop("scale parameter must be fixed for generalized binary priors.")
  }

  if (fix_g) {
    g <- g_init
    llik <- calc_gb_posterior(x, s, g, output = "llik")
  } else {
    s2 <- s^2

    ### separately optimize over parameters in different intervals for pi and
    ###   take the best solution
    opt_list <- list(NULL)
    val_list <- rep(NA, length(wlist) - 1)
    for (k in 1:(length(wlist) - 1)) {
      if (!is.null(g_init) &&
          g_init$pi[2] >= wlist[k] &&
          g_init$pi[2] <= wlist[k + 1]) {
        w <- g_init$pi[2]
        mu <- g_init$mean[2]
      } else {
        w <- (wlist[k] + wlist[k + 1]) / 2
        mu <- mu_init
      }

      ### update g
      iter <- 1
      continue_loop <- TRUE
      while (iter <= maxiter && continue_loop) {
        g <- tnormalmix(
          pi = c(1 - w, w),
          a = c(-Inf, 0),
          b = c(Inf, Inf),
          mean = c(0, mu),
          sd = c(0, abs(mu * scale))
        )
        zeta <- calc_gb_posterior(x, s, g, output = "zeta")

        # update g given the expectation of latent variables z
        w_new <- pmin(pmax(mean(zeta), wlist[k]), wlist[k+1])

        if (!identical(mode, "estimate")) {
          mu_new <- mode
        } else {
          # define the objective function to maximize to obtain mu
          opt_fn <- function(par) {
            par_sigma2 <- (par * scale)^2
            # calculate mean and variance of posterior corresponding to the
            #   truncated normal component in the mixture prior
            wgts <- s2 / (s2 + par_sigma2)
            post_mean <- wgts * par + (1 - wgts) * x
            post_var <- 1 / (1 / par_sigma2 + 1 / s2)
            post_sd <- pmax(sqrt(post_var), 1e-15)
            llik <- -0.5 * (log(par_sigma2 + s2) + (x - par)^2 / (par_sigma2 + s2)) +
              pnorm(-post_mean / post_sd, lower.tail = FALSE, log.p = TRUE) -
              pnorm(-par / sqrt(par_sigma2), lower.tail = FALSE, log.p = TRUE)
            # calculate objective function
            return(-sum(zeta * llik))
          }

          # update mu
          mu_new <- do.call(
            optimize,
            c(list(f = opt_fn,
                   interval = mu_range),
              control)
          )$minimum
        }

        # check for stopping criterion
        if (abs(w_new - w) < tol & abs(mu_new - mu) < mu * tol) {
          continue_loop <- FALSE
        }

        mu <- mu_new
        w <- w_new
        iter <- iter + 1
      }

      ### return the estimated g and likelihood
      g <- tnormalmix(
        pi = c(1 - w, w),
        a = c(-Inf, 0),
        b = c(Inf, Inf),
        mean = c(0, mu),
        sd = c(0, abs(mu * scale))
      )
      opt_list[[k]] <- g
      val_list[k] <- calc_gb_posterior(x, s, g, output = "llik")
    }

    ### take the optimal g among all intervals
    g <- opt_list[[which.max(val_list)]]
    llik <- max(val_list)
  }

  # Compute results.
  retlist <- list()

  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, x, s)
  }

  if (posterior_in_output(output)) {
    posterior <- calc_gb_posterior(x, s, g, output = "posterior")
    retlist <- add_posterior_to_retlist(retlist, posterior, output, x)
  }

  if (g_in_output(output)) {
    retlist <- add_g_to_retlist(retlist, g)
  }

  if (llik_in_output(output)) {
    retlist <- add_llik_to_retlist(retlist, llik, x, df = 2 * (1 - fix_g))
  }

  if (sampler_in_output(output)) {
    sampler <- calc_gb_posterior(x, s, g, output = "sampler")
    retlist <- add_sampler_to_retlist(retlist, sampler)
  }

  return(retlist)
}

#' @importFrom truncnorm etruncnorm vtruncnorm rtruncnorm
#' @importFrom stats rbinom
#'
calc_gb_posterior <- function(x, s, g, output){
  s2 <- s^2
  w <- g$pi[2]
  mu <- g$mean[2]
  sigma2 <- g$sd[2]^2

  # Calculate mean and variance of posterior corresponding to the truncated
  #   normal component in the mixture prior:
  wgts <- s2 / (s2 + sigma2)
  post_mean <- wgts * mu + (1 - wgts) * x
  post_var <- 1 / (1 / sigma2 + 1 / s2)
  post_sd <- pmax(sqrt(post_var), 1e-15)

  # Calculate marginal log likelihood for each mixture component:
  llik_mat <- matrix(NA, nrow = length(x), ncol = 2)
  llik_mat[,1] <- -0.5 * (log(s2) + x^2 / s2)
  llik_mat[,2] <- -0.5 * (log(sigma2 + s2) + (x - mu)^2 / (sigma2 + s2)) +
    pnorm(-post_mean / post_sd, lower.tail = FALSE, log.p = TRUE) -
    pnorm(-mu / sqrt(sigma2), lower.tail = FALSE, log.p = TRUE)
  llik_norms <- apply(llik_mat, 1, max)
  L_mat <- exp(llik_mat - llik_norms)

  # Calculate posterior weight for truncated normal in the posterior:
  zeta_mat <- t(t(L_mat) * g$pi)
  zeta_mat <- zeta_mat * (1 / rowSums(zeta_mat))
  zeta <- zeta_mat[,2]

  out <- list()
  if (output == "zeta") {
    out <- zeta
  } else if (output == "llik") {
    out <- sum(log(L_mat %*% g$pi)) + sum(llik_norms) - 0.5 * log(pi) * length(x)
  } else if (output == "posterior") {
    tmp1 <- etruncnorm(a = 0, mean = post_mean, sd = post_sd)
    tmp1[is.nan(tmp1)] <- 0
    tmp1[tmp1 < 0] <- 0

    tmp2 <- vtruncnorm(a = 0, mean = post_mean, sd = post_sd)
    tmp2[is.nan(tmp2)] <- 0
    tmp2[tmp2 < 0] <- 0

    pm <- zeta * tmp1
    pm2 <- zeta * (tmp1^2 + tmp2)
    pm2 <- pmax(pm2, pm^2)
    psd <- sqrt(pmax(0, pm2 - pm^2))
    lfsr <- 1 - zeta

    out <- list(mean = pm, mean2 = pm2, sd = psd, lfsr = lfsr)
  } else if (output == "sampler") {
    out <- function(nsamp) {
      nobs <- length(x)
      is_nonnull <- rbinom(nsamp * nobs, 1, rep(zeta, each = nsamp))
      samp <- is_nonnull * rtruncnorm(nsamp * nobs,
                                      a = 0,
                                      b = Inf,
                                      mean = rep(post_mean, each = nsamp),
                                      sd = rep(post_sd, each = nsamp))
      samp <- matrix(samp, nrow = nsamp)
      colnames(samp) <- names(x)
      return(samp)
    }
  }

  return(out)
}
