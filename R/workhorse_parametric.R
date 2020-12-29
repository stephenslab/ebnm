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
                                 use_grad,
                                 use_hess,
                                 checkg_fn,
                                 gtopar_fn,
                                 precomp_fn,
                                 nllik_fn,
                                 postcomp_fn,
                                 summres_fn,
                                 partog_fn,
                                 postsamp_fn,
                                 call) {
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
  par_init <- do.call(gtopar_fn, list(g_init = g_init,
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
                           startpar_fn = startpar_fn,
                           precomp_fn = precomp_fn,
                           nllik_fn = nllik_fn,
                           postcomp_fn = postcomp_fn,
                           optmethod = optmethod,
                           control = control,
                           use_grad = use_grad,
                           use_hess = use_hess)

  # Build return object.
  retlist <- list()
  if (posterior_in_output(output)) {
    posterior <- do.call(summres_fn, list(x = x,
                                          s = s,
                                          optpar = optres$par,
                                          output = output))
    retlist   <- add_posterior_to_retlist(retlist, posterior, output)
  }
  if (g_in_output(output)) {
    fitted_g <- do.call(partog_fn, list(par = optres$par))
    retlist  <- add_g_to_retlist(retlist, fitted_g)
  }
  if (llik_in_output(output)) {
    loglik  <- optres$val
    retlist <- add_llik_to_retlist(retlist, loglik)
  }
  if (sampler_in_output(output)) {
    post_sampler <- function(nsamp) {
      postsamp_fn(x, s, optres$par, nsamp)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

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
                         scale_name) {
  if (!is.null(g_init)) {
    if (!inherits(g_init, class_name)) {
      stop("g_init must be NULL or an object of class ", class_name, ".")
    }

    ncomp <- length(g_init$pi)
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
          && !isTRUE(all.equal(g_init$mean[1], mode))) {
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
