mle_parametric <- function(x,
                           s,
                           par,
                           fix_par,
                           precomp_fn,
                           nllik_fn,
                           processoptres_fn,
                           optmethod,
                           control,
                           use_grad,
                           use_hess) {
  precomp  <- do.call(precomp_fn, list(x = x, s = s, par = par, fix_par = fix_par))

  # Parameters that end up getting passed to all optimization functions.
  fn_params <- c(list(x = x, s = s, par = par, fix_par = fix_par), precomp)

  if (optmethod == "nlm") {
    control <- modifyList(nlm_control_defaults(), control)

    optres <- do.call(nlm, c(list(f = nllik_fn, p = startpar), fn_params,
                             list(calc_grad = use_grad, calc_hess = use_hess),
                             control))
    opttheta <- optres$estimate
    optval   <- optres$minimum
  } else if (optmethod == "trust") {
    control <- modifyList(trust_control_defaults(), control)

    # trust requires a gradient and Hessian.
    fn <- function(theta, ...) {
      nllik <- do.call(nllik_fn,
                       list(theta = theta, calc_grad = TRUE, calc_hess = TRUE, ...))
      return(list(value = nllik,
                  gradient = attr(nllik, "gradient"),
                  hessian = attr(nllik, "hessian")))
    }

    optres <- do.call(trust, c(list(objfun = fn, parinit = startpar), fn_params,
                               control))
    opttheta <- optres$argument
    optval   <- optres$value
  } else if (optmethod == "lbfgsb") {
    control <- modifyList(lbfgsb_control_defaults(), control)

    # optim cannot accept a Hessian.
    fn <- function(theta, ...) {
      return(do.call(nllik_fn,
                     list(theta = theta, calc_grad = FALSE, calc_hess = FALSE, ...)))
    }
    if (use_grad) {
      gr <- function(theta, ...) {
        nllik <- do.call(nllik_fn,
                         list(theta = theta, calc_grad = TRUE, calc_hess = FALSE, ...))
        return(attr(nllik, "gradient"))
      }
    } else {
      gr <- NULL
    }

    optres  <- do.call(optim, c(list(par = startpar, fn = fn, gr = gr), fn_params,
                                control, list(method = "L-BFGS-B")))
    opttheta  <- optres$par
    optval    <- optres$value
  } else if (optmethod == "optimize") {
    control <- modifyList(optimize_control_defaults(), control)

    optres <- do.call(optimize, c(list(f = nllik_fn), fn_params,
                                  list(calc_grad = FALSE, calc_hess = FALSE),
                                  control))
    opttheta <- optres$minimum
    optval   <- optres$objective
  }

  retlist <- do.call(processoptres_fn, c(list(opttheta = opttheta,
                                              optval = optval,
                                              par = par,
                                              fix_par = fix_par),
                                         precomp))

  return(retlist)
}
