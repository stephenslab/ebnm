#' Solve the EBNM problem for grouped data
#'
#' Solves the empirical Bayes normal means (EBNM) problem for observations
#'   belonging to distinct groups.
#'
#' The EBNM model for grouped data, with observations \eqn{x_j} belonging to
#'   groups \eqn{k = 1, ..., K}, is
#'   \deqn{x_j | \theta_j, s_j \sim N(\theta_j, s_j^2)}
#'   \deqn{\theta_j \sim g_{k(j)} \in G_{k(j)}.}
#'
#' Solving the EBNM problem for grouped data is equivalent to solving a
#'   separate EBNM problem for each group \eqn{k = 1, ..., K}, with the optimal
#'   log likelihood equal to the sum of the optimal log likelihoods for each
#'   separate problem.
#'
#' @inheritParams ebnm
#' @inherit ebnm return
#'
#' @param group A vector of character strings that gives the group to which each
#'   observation belongs. It must have the same length as argument \code{x}. For
#'   an example of usage, see Examples below.
#'
#' @param prior_family A named vector that specifies the prior family \eqn{G}
#'   for each group. If the same prior family is to be used for all groups, then
#'   a character string may be used instead.
#'
#' @param mode A named list that specifies, for each group, the mode of the
#'   respective prior \eqn{g}, or \code{"estimate"} if the mode is to be
#'   estimated from the data. If the mode is the same across groups, then a
#'   scalar may be used instead. If all modes are to be estimated, then
#'   \code{mode = "estimate"} may be used.
#'
#' @param scale A named list that specifies, for each group, the scale
#'   parameter(s) of the respective prior, or \code{"estimate"} if the scale
#'   parameters are to be estimated from the data. If the scale parameter is the
#'   same across groups, then a scalar may be used instead. If all scales are to
#'   be estimated, then \code{scale = "estimate"} may be used.
#'
#' @param g_init The prior distributions \eqn{g}. Usually this is left
#'   unspecified (\code{NULL}) and estimated from the data. However, it can be
#'   used in conjuction with \code{fix_g = TRUE} to fix the prior (useful, for
#'   example, to do computations with the "true" \eqn{g} in simulations). If
#'   \code{g_init} is specified but \code{fix_g = FALSE}, \code{g_init}
#'   specifies the initial value of \eqn{g} used during optimization. If
#'   \code{g_init} is supplied, it should be a named list that specifies, for
#'   each group, a prior of the appropriate class (\code{\link[ashr]{normalmix}}
#'   for normal, point-normal,
#'   scale mixture of normals, and \code{deconvolveR} prior families, as well as
#'   for the NPMLE; class \code{\link{laplacemix}} for
#'   point-Laplace families; class \code{\link{gammamix}} for point-exponential
#'   families; class \code{\link{horseshoe}} for horseshoe families; and class
#'   \code{\link[ashr]{unimix}} for \code{unimodal_} families).
#'
#' @seealso \code{\link{ebnm}}
#'
#' @examples
#' group <- c(rep("small_sd", 100), rep("large_sd", 100))
#' theta <- c(rnorm(100, sd = 1), rnorm(100, sd = 10))
#' s <- 1
#' x <- theta + rnorm(200, 0, s)
#'
#' ebnm.group.res <- ebnm_group(x, s, group)
#'
#' # Use different prior families for each group:
#' ebnm.group.res <- ebnm_group(
#'   x, s, group,
#'   prior_family = list(small_sd = "normal", large_sd = "normal_scale_mixture")
#' )
#'
#' # Different modes and scales can be set similarly:
#' ebnm.group.res <- ebnm_group(
#'   x, s, group,
#'   mode = list(small_sd = 0, large_sd = "estimate"),
#'   scale = list(small_sd = 1, large_sd = "estimate")
#' )
#'
#' @export
#'
ebnm_group <- function(x,
                       s = 1,
                       group,
                       prior_family = "point_normal",
                       mode = 0,
                       scale = "estimate",
                       g_init = NULL,
                       fix_g = FALSE,
                       output = ebnm_output_default(),
                       ...) {
  check_args_group(x, group, prior_family, mode, scale, g_init)

  args_list <- sapply(unique(group), function(grp) {
    idx = which(group == grp)
    grp_args <- list(
      group = grp,
      x = x[idx],
      s = ifelse(length(s) == 1, s, s[idx]),
      prior_family = ifelse(is.list(prior_family), prior_family[[grp]], prior_family),
      mode = ifelse(is.list(mode), mode[[grp]], mode),
      scale = ifelse(is.list(scale), scale[[grp]], scale),
      fix_g = fix_g,
      output = output
    )
    if (!is.null(g_init)) {
      grp_args[["g_init"]] <- g_init[[grp]]
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
    lliks   <- lapply(ebnm_res, `[[`, llik_ret_str())
    loglik  <- as.numeric(Reduce(`+`, lliks))
    df      <- sum(sapply(lliks, attr, "df"))
    retlist <- add_llik_to_retlist(retlist, loglik, x, df)
  }

  # TODO: There is as yet no way to set parameter burn for ebnm_horseshoe.
  if (sampler_in_output(output)) {
    post_sampler <- function(nsamp) {
      samp <- matrix(nrow = nsamp, ncol = length(x))
      samp_res <- sapply(
        ebnm_res,
        function(grp) grp[[samp_ret_str()]](nsamp),
        simplify = FALSE
      )
      for (grp in unique(group)) {
        idx <- which(group == grp)
        samp[, idx] <- samp_res[[grp]]
      }
      colnames(samp) <- names(x)
      return(samp)
    }
    retlist <- add_sampler_to_retlist(retlist, post_sampler)
  }

  return(as_ebnm(retlist, match.call()))
}

check_args_group <- function(x, group, prior_family, mode, scale, g_init) {
  if (!(length(group) == length(x))) {
    stop("Argument 'group' must be the same length as argument 'x'.")
  }

  if (any(is.na(group))) {
    stop("Missing group assignments are not allowed.")
  }

  grps <- unique(group)

  check_named_list_arg(prior_family, grps, "prior_family")
  check_named_list_arg(mode, grps, "mode")
  check_named_list_arg(scale, grps, "scale")

  if (!is.null(g_init)) {
    if (!(is.list(g_init) && all(grps %in% names(g_init)))) {
      stop("Argument 'g_init' must either be NULL or a named list with names ",
           "corresponding to the unique values in argument 'group'.")
    }
  }
}

check_named_list_arg <- function(arg, grps, arg_name) {
  if (!(is.list(arg) && all(grps %in% names(arg)))) {
    if (length(arg) != 1) {
      stop("Argument ", arg_name, " must either have length 1 or be a named list ",
           "with names corresponding to the unique values in argument 'group'.")
    }
  }
}
