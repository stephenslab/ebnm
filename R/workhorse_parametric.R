# startpar_fn provides initial values for parameters that are to be estimated
#   and takes arguments x, s, par, and fix_par.
#
# precomp_fn does precomputations (if any). It also takes arguments x, s, par,
#   and fix_par.
#
# nllik_fn calculates the negative log likelihood, gradient (if
#   calc_grad = TRUE), and Hessian (if calc_hess = TRUE). It takes arguments
#   x, s, par, fix_par, calc_grad, and calc_hess, and then whatever has been
#   calculated by precomp_fn.
#
# gfromopt_fn converts the optimization results into a suitable return object
#   (transforming parameters when necessary). It takes arguments opttheta, g,
#   and fix_par, and then whatever has been calculated by precomp_fn.

#
parametric_workhorse <- function(x,
                                 s,
                                 mode,
                                 scale,
                                 pointmass,
                                 g_init,
                                 fix_g,
                                 output,
                                 optmethod,
                                 control,
                                 checkg_fn,
                                 initpar_fn,
                                 scalepar_fn,
                                 precomp_fn,
                                 nllik_fn,
                                 postcomp_fn,
                                 summres_fn,
                                 partog_fn,
                                 postsamp_fn,
                                 call) {
  # I'm not sure why this is the case, but I run into infinite recursion issues
  #   if I don't extract the things I need from call here.
  call <- list(mode = call$mode, scale = call$scale)

  # Check that argument g_init is valid. All parametric families currently
  #   call into function check_g_init below.
  do.call(checkg_fn, list(g_init = g_init,
                          fix_g = fix_g,
                          mode = mode,
                          scale = scale,
                          pointmass = pointmass,
                          call = call))

  # Translate ebnm interface (mode/scale/pointmass/g_init/fix_g) into a
  #   generalized optimization interface (par_init/fix_par). fix_par, a vector
  #   of length 3, indicates whether 1. the weight of the spike component;
  #   2. the scale of the slab component; and 3. the location of the components
  #   is fixed.
  par_init <- do.call(initpar_fn, list(g_init = g_init,
                                       mode = mode,
                                       scale = scale,
                                       pointmass = pointmass,
                                       x = x,
                                       s = s))
  if (fix_g) {
    fix_par <- c(TRUE, TRUE, TRUE)
  } else {
    fix_par <- c(!pointmass,
                 !identical(scale, "estimate"),
                 !identical(mode, "estimate"))
  }

  optmethod <- handle_optmethod_parameter(optmethod, fix_par)

  # Don't use observations with infinite SEs when estimating g.
  x_optset <- x
  s_optset <- s
  if (any(is.infinite(s))) {
    x_optset <- x[is.finite(s)]
    s_optset <- s[is.finite(s)]
  }

  # Estimate g. Function mle_parametric returns a list with fields par (which
  #   gives the estimated values of the parameters) and val (the optimal
  #   log likelihood attained).
  optres <- mle_parametric(x = x_optset,
                           s = s_optset,
                           par_init = par_init,
                           fix_par = fix_par,
                           scalepar_fn = scalepar_fn,
                           precomp_fn = precomp_fn,
                           nllik_fn = nllik_fn,
                           postcomp_fn = postcomp_fn,
                           optmethod = optmethod$fn,
                           control = control,
                           use_grad = optmethod$use_grad,
                           use_hess = optmethod$use_hess)

  # Build return object.
  retlist <- list()

  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, x, s)
  }

  if (posterior_in_output(output)) {
    posterior <- do.call(summres_fn, list(x = x,
                                          s = s,
                                          optpar = optres$par,
                                          output = output))
    retlist <- add_posterior_to_retlist(retlist, posterior, output, x)
  }

  if (g_in_output(output)) {
    fitted_g <- do.call(partog_fn, list(par = optres$par))
    retlist  <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    loglik  <- optres$val
    retlist <- add_llik_to_retlist(retlist, loglik, x, df = sum(!fix_par))
  }

  if (sampler_in_output(output)) {
    sampler <- function(nsamp) {
      samp <- postsamp_fn(x, s, optres$par, nsamp)
      colnames(samp) <- names(x)
      return(samp)
    }
    retlist <- add_sampler_to_retlist(retlist, sampler)
  }

  return(retlist)
}


handle_optmethod_parameter <- function(optmethod, fix_par) {
  optmethod <- match.arg(optmethod, c("nohess_nlm", "nlm", "nograd_nlm",
                                      "lbfgsb", "nograd_lbfgsb",
                                      "trust",
                                      "optimize"))
  return(
    switch(optmethod,
           nlm = list(fn = "nlm", use_grad = TRUE, use_hess = TRUE),
           lbfgsb = list(fn = "lbfgsb", use_grad = TRUE, use_hess = FALSE),
           trust = list(fn = "trust", use_grad = TRUE, use_hess = TRUE),
           nograd_nlm = list(fn = "nlm", use_grad = FALSE, use_hess = FALSE),
           nograd_lbfgsb = list(fn = "lbfgsb", use_grad = FALSE, use_hess = FALSE),
           nohess_nlm = list(fn = "nlm", use_grad = TRUE, use_hess = FALSE),
           optimize = list(fn = "optimize", use_grad = FALSE, use_hess = FALSE))
  )
}


# This function calls the selected optimization routine.
#
#' @importFrom stats median nlm optim optimize
#' @importFrom trust trust
#' @importFrom utils modifyList
#'
mle_parametric <- function(x,
                           s,
                           par_init,
                           fix_par,
                           scalepar_fn,
                           precomp_fn,
                           nllik_fn,
                           postcomp_fn,
                           optmethod,
                           control,
                           use_grad,
                           use_hess) {
  scale_factor <- 1 / median(s[s > 0])
  x <- x * scale_factor
  s <- s * scale_factor

  par_init <- do.call(scalepar_fn, list(par = par_init,
                                        scale_factor = scale_factor))

  precomp <- do.call(precomp_fn, list(x = x,
                                      s = s,
                                      par_init = par_init,
                                      fix_par = fix_par))

  # Parameters that end up getting passed to all optimization functions:
  fn_params <- c(list(x = x, s = s, par_init = par_init, fix_par = fix_par),
                 precomp)

  p <- unlist(par_init)[!fix_par]

  # Fix issue #46 (don't initialize using pi0 = 0 or pi0 = 1):
  if (!(fix_par[1]) && is.infinite(p[1])) {
    p[1] <- sign(p[1]) * log(length(x))
  }

  if (all(fix_par)) {
    optpar <- par_init
    optval <- do.call(nllik_fn, c(list(par = NULL), fn_params,
                                  list(calc_grad = FALSE, calc_hess = FALSE)))
  } else if (optmethod == "nlm") {
    control <- modifyList(nlm_control_defaults(), control)

    optres <- do.call(nlm, c(list(f = nllik_fn, p = p),
                             fn_params,
                             list(calc_grad = use_grad, calc_hess = use_hess),
                             control))
    optpar <- optres$estimate
    optval <- optres$minimum
  } else if (optmethod == "trust") {
    control <- modifyList(trust_control_defaults(), control)

    # trust requires both a gradient and a Hessian.
    fn <- function(par, ...) {
      nllik <- do.call(nllik_fn,
                       list(par = par, calc_grad = TRUE, calc_hess = TRUE, ...))
      return(list(value = nllik,
                  gradient = attr(nllik, "gradient"),
                  hessian = attr(nllik, "hessian")))
    }
    optres <- do.call(trust::trust, c(list(objfun = fn, parinit = p),
                               fn_params,
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

    optres <- do.call(optim, c(list(par = p, fn = fn, gr = gr),
                               fn_params,
                               list(control = control),
                               list(method = "L-BFGS-B")))
    optpar <- optres$par
    optval <- optres$value
  } else if (optmethod == "optimize") {
    control <- modifyList(optimize_control_defaults(), control)

    optres <- do.call(optimize, c(list(f = nllik_fn), fn_params,
                                  list(calc_grad = FALSE, calc_hess = FALSE),
                                  control))
    optpar <- optres$minimum
    optval <- optres$objective
  }

  # Combine the fixed and estimated parameters.
  retpar <- par_init
  retpar[!fix_par] <- optpar

  # Re-scale parameters and log likelihood.
  retpar <- do.call(scalepar_fn, list(par = retpar,
                                      scale_factor = 1 / scale_factor))
  optval <- optval - sum(is.finite(x)) * log(scale_factor)

  retlist <- do.call(postcomp_fn, c(list(optpar = retpar,
                                         optval = optval,
                                         x = x,
                                         s = s,
                                         par_init = par_init,
                                         fix_par = fix_par,
                                         scale_factor = scale_factor),
                                    precomp))

  return(retlist)
}


# This function is currently used by all parametric families to check g_init.
#
check_g_init <- function(g_init,
                         fix_g,
                         mode,
                         scale,
                         pointmass,
                         call,
                         class_name,
                         scale_name,
                         mode_name = "mean") {
  if (!is.null(g_init)) {
    if (!inherits(g_init, class_name)) {
      stop("g_init must be NULL or an object of class ", class_name, ".")
    }

    ncomp <- length(g_init[[1]])
    if (!(ncomp == 1 || (pointmass && ncomp == 2))) {
      stop("g_init does not have the correct number of components.")
    }
    if (ncomp == 2 && g_init[[scale_name]][1] != 0) {
      stop("The first component of g_init must be a point mass.")
    }

    if (fix_g && (!is.null(call$mode) || !is.null(call$scale))) {
      warning("mode and scale parameters are ignored when g is fixed.")
    }

    if (!fix_g) {
      # all.equal allows for numerical error:
      if (!is.null(call$mode)
          && !identical(mode, "estimate")
          && !isTRUE(all.equal(g_init[[mode_name]], rep(mode, length(g_init[[mode_name]]))))) {
        stop("If mode is fixed and g_init is supplied, they must agree.")
      }
      g_scale <- g_init[[scale_name]][ncomp]
      if (!is.null(call$scale)
          && !identical(scale, "estimate")
          && !isTRUE(all.equal(g_scale, scale))) {
        stop("If scale is fixed and g_init is supplied, they must agree.")
      }
    }
  }
}
