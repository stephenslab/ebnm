#' @importFrom stats nlm
#'
mle_point_laplace <- function(x, s, g, control,
                              optmethod, use_grad, use_hess) {
  startpar <- pl_startpar(x, s, g)

  lf <- calc_lf(x, s)

  fn_params <- list(x = x, s = s, lf = lf)

  if (optmethod == "nlm") {
    control <- modifyList(nlm_control_defaults(), control)

    optres <- do.call(nlm, c(list(f = pl_nllik, p = startpar), fn_params,
                             list(calc_grad = use_grad, calc_hess = use_hess),
                             control))
    optpar <- optres$estimate
    optval <- optres$minimum
  } else if (optmethod == "trust") {
    control <- modifyList(trust_control_defaults(), control)

    # trust requires a gradient and Hessian.
    fn <- function(par, ...) {
      nllik <- pl_nllik(par, calc_grad = TRUE, calc_hess = TRUE, ...)
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
      return(pl_nllik(par, calc_grad = FALSE, calc_hess = FALSE, ...))
    }
    if (use_grad) {
      gr <- function(par, ...) {
        nllik <- pl_nllik(par, calc_grad = TRUE, calc_hess = FALSE, ...)
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

  retlist <- pl_g_from_optpar(optpar)
  retlist$val <- -optval

  # Check the solution pi0 = 1.
  if (sum(lf) > retlist$val) {
    retlist$pi0 <- 1
    retlist$a <- 1
    retlist$val <- sum(lf)
  }

  return(retlist)
}

# Initial values.
pl_startpar <- function(x, s, g) {
  startpar <- numeric(0)

  if (!is.null(g$pi0) && g$pi0 > 0 && g$pi0 < 1) {
    startpar <- c(startpar, log(1 / g$pi0 - 1))
  } else {
    startpar <- c(startpar, 0) # default for -logit(pi0)
  }

  if (!is.null(g$a)) {
    startpar <- c(startpar, log(g$a))
  } else {
    startpar <- c(startpar, -0.5 * log(mean(x^2) / 2)) # default for log(a)
  }

  return(startpar)
}

# Pull pi0 and a out of the optimization results.
pl_g_from_optpar <- function(par) {
  return(list(pi0 = 1 / (1 + exp(par[1])), a = exp(par[2])))
}
