#' @importFrom stats nlm
#'
mle_point_normal <- function(x, s, g, control, fix_pi0, fix_a, fix_mu) {
  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }

  startpar <- pn_startpar(x, s, g, fix_pi0, fix_a, fix_mu)

  if (fix_pi0) {
    alpha <- -log(1 / g$pi0 - 1)
  } else {
    alpha <- NULL
  }
  if (fix_a) {
    beta <- -log(g$a)
  } else {
    beta <- NULL
  }
  if (fix_mu) {
    mu <- g$mu
  } else {
    mu <- NULL
  }

  if (any(s == 0)) {
    which_s0 <- which(s == 0)
    which_x_nz <- which(x[which_s0] != mu)
    n0 <- length(which_s0) - length(which_x_nz)
    n1 <- length(which_x_nz)
    sum1 <- sum((x[which_s0[which_x_nz]] - mu)^2)
    x <- x[-which_s0]
    s <- s[-which_s0]
  } else {
    n0 <- 0
    n1 <- 0
    sum1 <- 0
  }
  n2 <- length(x)

  s2 <- s^2

  if (fix_mu) {
    z <- (x - mu)^2 / s2
    sum_z <- sum(z)
  } else {
    z <- NULL
    sum_z <- NULL
  }

  fn_params <- list(fix_pi0 = fix_pi0, fix_a = fix_a, fix_mu = fix_mu,
                    alpha = alpha, beta = beta, mu = mu,
                    n0 = n0, n1 = n1, sum1 = sum1, n2 = n2,
                    x = x, s2 = s2, z = z, sum_z = sum_z)

  # Sometimes nlm thinks that the gradient is being calculated incorrectly.
  #   Reducing the number of significant digits often solves the problem.
  control <- modifyList(list(ndigit = 8), control)

  optres <- do.call(nlm, c(list(pn_nlm_fn, startpar), fn_params, control))

  # if (inherits(optres, "try-error")) {
  #   stop("Re-attempt at optimization failed, possibly due to one or more ",
  #        "very small standard errors. The smallest nonzero standard error ",
  #        "seen during optimization was ", signif(min(s), 2), ".")
  # }

  retlist <- pn_g_from_optpar(optres$estimate, g, fix_pi0, fix_a, fix_mu)
  retlist$val <- pn_llik_from_optval(optres$minimum, n1, n2, s2)

  # Check the solution pi0 = 1.
  if (!fix_pi0 && fix_mu && n1 == 0) {
    pi0_val <- pn_llik_from_optval(sum_z / 2, n1, n2, s2)
    if (pi0_val > retlist$val) {
      retlist$pi0 <- 1
      retlist$a <- 1
      retlist$val <- pi0_val
    }
  }

  return(retlist)
}

# Initial values.
pn_startpar <- function(x, s, g, fix_pi0, fix_a, fix_mu) {
  startpar <- numeric(0)

  if (!fix_pi0) {
    if (!is.null(g$pi0) && g$pi0 > 0 && g$pi0 < 1) {
      startpar <- c(startpar, -log(1 / g$pi0 - 1))
    } else {
      startpar <- c(startpar, 0) # default for logit(pi0)
    }
  }

  if (!fix_a) {
    if (!is.null(g$a)) {
      startpar <- c(startpar, -log(g$a))
    } else {
      startpar <- c(startpar, log(mean(x^2))) # default for -log(a)
    }
  }

  if (!fix_mu) {
    if (!is.null(g$mu)) {
      startpar <- c(startpar, g$mu)
    } else {
      startpar <- c(startpar, mean(x)) # default for mu
    }
  }

  return(startpar)
}

# Pull pi0, a, and mu out of the optimization results.
pn_g_from_optpar <- function(optpar, g, fix_pi0, fix_a, fix_mu) {
  opt_g <- list()

  i <- 1
  if (fix_pi0) {
    opt_g$pi0 <- g$pi0
  } else {
    opt_g$pi0 <- 1 / (exp(-optpar[i]) + 1)
    i <- i + 1
  }

  if (fix_a) {
    opt_g$a <- g$a
  } else {
    opt_g$a <- exp(-optpar[i])
    i <- i + 1
  }

  if (fix_mu) {
    opt_g$mu <- g$mu
  } else {
    opt_g$mu <- optpar[i]
  }

  return(opt_g)
}

pn_llik_from_optval <- function(optval, n1, n2, s2) {
  if (length(s2) == 1) {
    sum.log.s2 <- n2 * log(s2)
  } else {
    sum.log.s2 <- sum(log(s2))
  }
  return(-optval - 0.5 * ((n1 + n2) * log(2 * pi) + sum.log.s2))
}
