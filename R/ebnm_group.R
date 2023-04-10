ebnm_group <- function(x,
                       s = 1,
                       group,
                       prior_family = "point_normal",
                       mode = 0,
                       scale = "estimate",
                       g_init = NULL,
                       fix_g = FALSE,
                       output = output_default(),
                       ...) {
  # Check group parameter.
  # Check prior_family parameter.
  # Check g_init.

  args_list <- sapply(unique(group), function(grp) {
    idx = which(group == grp)
    grp_args <- list(
      x = x[idx],
      s = ifelse(length(s) == 1, s, s[idx]),
      prior_family = ifelse(length(prior_family) == 1, prior_family, prior_family[grp]),
      g_init = NULL,
      fix_g = fix_g,
      output = output
    )
    if (!is.null(g_init)) {
      grp_args$g_init <- g_init[[grp]]
    }
    return(grp_args)
  }, simplify = FALSE)

  ebnm_res <- sapply(args_list, function(args) {
    grp <- args$group
    args$group <- NULL
    tryCatch(
      res <- do.call(ebnm, args),
      error = function(e) stop(paste0("Error in group ", grp, ": ", e$message))
    )
    return(res)
  }, simplify = FALSE)

  # Compute results.
  retlist <- list()

  if (data_in_output(output)) {
    retlist <- add_data_to_retlist(retlist, x, rep(s, length.out = length(x)))
    retlist[[data_ret_str()]][[grp_ret_str()]] <- group
  }

  if (posterior_in_output(output)) {
    fields <- names(ebnm_res[[1]][[df_ret_str()]])

    # Initialize empty data frame.
    posterior <- data.frame(sapply(fields, function(field) rep(0, length(x))))

    # Loop over groups and fields.
    for (grp in unique(group)) {
      idx <- which(group == grp)
      for (field in fields) {
        posterior[[field]][idx] <- ebnm_res[[grp]][[df_ret_str()]][[field]]
      }
    }

    retlist[[df_ret_str()]] <- posterior
  }

  if (g_in_output(output)) {
    fitted_g <- sapply(ebnm_res, `[[`, g_ret_str(), simplify = FALSE)
    retlist <- add_g_to_retlist(retlist, fitted_g)
  }

  if (llik_in_output(output)) {
    loglik  <- sum(sapply(ebnm_res, `[[`, llik_ret_str()))
    retlist <- add_llik_to_retlist(retlist, loglik)
  }

  # TODO: Posterior sampler. Handle ebnm_horseshoe.
  if (sampler_in_output(output)) {
    post_sampler <- function(nsamp) {
      samp <- matrix(nrow = nsamp, ncol = length(x))
      samp_res <- sapply(
        ebnm_res,
        function (grp) grp[[samp_ret_str()]](nsamp),
        simplify = FALSE
      )
      for (grp in unique(group)) {
        idx <- which(group == grp)
        samp[, idx] <- samp_res[[grp]]
      }
      return(samp)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(as_ebnm(retlist, match.call()))
}
