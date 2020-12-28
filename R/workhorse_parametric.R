parametric_workhorse <- function(x,
                                 s,
                                 pointmass,
                                 mode,
                                 scale,
                                 g_init,
                                 fix_g,
                                 output,
                                 optmethod,
                                 use_grad,
                                 use_hess,
                                 control,
                                 checkg_fn,
                                 gtopar_fn,
                                 startpar_fn,
                                 precomp_fn,
                                 nllik_fn,
                                 gfromopt_fn,
                                 summres_fn,
                                 partog_fn,
                                 postsamp_fn,
                                 call) {
  if (length(scale) != 1) {
    stop("Argument 'scale' must be either 'estimate' or a scalar.")
  }

  do.call(checkg_fn, list(g_init = g_init, fix_g = fix_g,
                          pointmass = pointmass, mode = mode, scale = scale,
                          call = call))

  # Assume three parameters (for now), which respectively parametrize the
  #   weight of the spike component, the scale of the slab component, and the
  #   location of the spike component.
  fix_par <- c(!pointmass,
               !identical(scale, "estimate"),
               !identical(mode, "estimate"))
  if (all(fix_par)) {
    fix_g <- TRUE
  }

  par <- do.call(gtopar_fn, list(g_init = g_init, fix_par = fix_par))

  x_optset <- x
  s_optset <- s

  # Don't use observations with infinite SEs when estimating g.
  if (any(is.infinite(s))) {
    x_optset <- x[is.finite(s)]
    s_optset <- s[is.finite(s)]
  }

  # Estimate g.
  optpar <- mle_parametric(x, s, par, fix_par,
                           startpar_fn, precomp_fn, nllik_fn, gfromopt_fn,
                           optmethod, control, use_grad, use_hess)

  retlist <- list()

  if (posterior_in_output(output)) {
    posterior <- do.call(summres_fn, list(x = x, s = s, g = optpar, output = output))
    retlist   <- add_posterior_to_retlist(retlist, posterior, output)
  }

  if (g_in_output(output)) {
    fitted_g <- do.call(partog_fn, list(par = optpar))
    retlist  <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    loglik <- optpar$val
    retlist <- add_llik_to_retlist(retlist, loglik)
  }

  if (sampler_in_output(output)) {
    post_sampler <- function(nsamp) {
      postsamp_fn(x, s, optpar, nsamp)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(retlist)
}

# Used by both ebnm_pn_workhorse and ebnm_pl_workhorse.
#
check_g_init <- function(g_init,
                         fix_g,
                         pointmass,
                         mode,
                         scale,
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
