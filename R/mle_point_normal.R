#' @importFrom stats nlm
#' @importFrom trust trust
#'
mle_point_normal <- function(x, s, g, control,
                             fix_pi0, fix_a, fix_mu,
                             optmethod, use_grad, use_hess) {
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

  if (optmethod == "nlm") {
    control <- modifyList(nlm_control_defaults(), control)

    optres <- do.call(nlm, c(list(f = pn_nllik, p = startpar), fn_params,
                             list(calc_grad = use_grad, calc_hess = use_hess),
                             control))
    optpar <- optres$estimate
    optval <- optres$minimum
  } else if (optmethod == "trust") {
    control <- modifyList(trust_control_defaults(), control)

    # trust requires a gradient and Hessian.
    fn <- function(par, ...) {
      nllik <- pn_nllik(par, calc_grad = TRUE, calc_hess = TRUE, ...)
      return(list(value = nllik,
                  gradient = attr(nllik, "gradient"),
                  hessian = attr(nllik, "hessian")))
    }

    optres <- do.call(trust, c(list(objfun = fn, parinit = startpar), fn_params,
                               control))
    optpar <- optres$argument
    optval <- optres$value
  } else if (optmethod == "lbfgsb") {
    control <- modifyList(lbfgsb_control_defaults(), control)

    # optim cannot accept a Hessian.
    fn <- function(par, ...) {
      return(pn_nllik(par, calc_grad = FALSE, calc_hess = FALSE, ...))
    }
    if (use_grad) {
      gr <- function(par, ...) {
        nllik <- pn_nllik(par, calc_grad = TRUE, calc_hess = FALSE, ...)
        return(attr(nllik, "gradient"))
      }
    } else {
      gr <- NULL
    }

    optres  <- do.call(optim, c(list(par = startpar, fn = fn, gr = gr), fn_params,
                                control, list(method = "L-BFGS-B")))
    optpar  <- optres$par
    optval  <- optres$value
  }

  retlist <- pn_g_from_optpar(optpar, g, fix_pi0, fix_a, fix_mu)
  retlist$val <- pn_llik_from_optval(optval, n1, n2, s2)

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
