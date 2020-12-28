# startpar_fn provides initial values for parameters that are to be estimated
#   and takes arguments x, s, g, and fix_par.
#
# precomp_fn does precomputations (if any). It also takes arguments x, s, g,
#   and fix_par.
#
# nllik_fn calculates the negative log likelihood, gradient (if
#   calc_grad = TRUE), and Hessian (if calc_hess = TRUE). It takes arguments
#   x, s, g, fix_par, calc_grad, and calc_hess, and then whatever has been
#   calculated by precomp_fn.
#
# gfromopt_fn converts the optimization results into a suitable return object
#   (transforming parameters when necessary). It takes arguments optpar, g,
#   and fix_par, and then whatever has been calculated by precomp_fn.

mle_parametric <- function(x, s, g, fix_par,
                           startpar_fn, precomp_fn, nllik_fn, gfromopt_fn,
                           optmethod, control, use_grad, use_hess) {
  startpar <- do.call(startpar_fn, list(x = x, s = s, g = g, fix_par = fix_par))
  precomp  <- do.call(precomp_fn, list(x = x, s = s, g = g, fix_par = fix_par))

  # Parameters that end up getting passed to all optimization functions.
  fn_params <- c(list(x = x, s = s, g = g, fix_par = fix_par), precomp)

  if (optmethod == "nlm") {
    control <- modifyList(nlm_control_defaults(), control)

    optres <- do.call(nlm, c(list(f = nllik_fn, p = startpar), fn_params,
                             list(calc_grad = use_grad, calc_hess = use_hess),
                             control))
    optpar <- optres$estimate
    optval <- optres$minimum
  } else if (optmethod == "trust") {
    control <- modifyList(trust_control_defaults(), control)

    # trust requires a gradient and Hessian.
    fn <- function(par, ...) {
      nllik <- do.call(nllik_fn,
                       list(par = par, calc_grad = TRUE, calc_hess = TRUE, ...))
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
      return(do.call(nllik_fn,
                     list(par = par, calc_grad = FALSE, calc_hess = FALSE, ...)))
    }
    if (use_grad) {
      gr <- function(par, ...) {
        nllik <- do.call(nllik_fn,
                         list(par = par, calc_grad = TRUE, calc_hess = FALSE, ...))
        return(attr(nllik, "gradient"))
      }
    } else {
      gr <- NULL
    }

    optres  <- do.call(optim, c(list(par = startpar, fn = fn, gr = gr), fn_params,
                                control, list(method = "L-BFGS-B")))
    optpar  <- optres$par
    optval  <- optres$value
  } else if (optmethod == "optimize") {
    control <- modifyList(optimize_control_defaults(), control)

    optres <- do.call(optimize, c(list(f = nllik_fn), fn_params,
                                  list(calc_grad = FALSE, calc_hess = FALSE),
                                  control))
    optpar <- optres$minimum
    optval <- optres$objective
  }

  retlist <- do.call(gfromopt_fn, c(list(optpar = optpar, optval = optval,
                                         g = g, fix_par = fix_par),
                                    precomp))

  return(retlist)
}
