#' Plot an ebnm object
#'
#' Given one or more fitted \code{\link{ebnm}} object(s), produces a plot of
#'   posterior means vs. observations. If desired, a plot of cumulative
#'   distribution functions of fitted prior(s) can also be produced.
#'
#' @param x The fitted \code{ebnm} object.
#'
#' @param incl_pm Plot posterior means vs. observations?
#'
#' @param incl_cdf Plot the cumulative distribution functions?
#'
#' @param subset The subset of observations to include on the plot of posterior
#'   means vs. observations. Can be a numeric vector corresponding to indices
#'   of observations to plot, or a character vector if observations are named.
#'   If \code{subset = NULL} then all observations will be plotted.
#'
#' @param remove_abline To better illustrate shrinkage effects, the plot of
#'   posterior means vs. observations includes the line \eqn{y = x} by default.
#'   If \code{remove_abline = TRUE}, then this line will not be drawn.
#'
#' @param ... Additional \code{ebnm} objects to be included on the same plots.
#'
#' @method plot ebnm
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @importFrom ggplot2 theme_minimal scale_color_brewer
#' @importFrom ggplot2 geom_line
#' @importFrom ashr mixcdf
#' @importFrom grDevices devAskNewPage
#'
#' @examples
#' theta <- c(rep(0, 100), rexp(100))
#' theta[1:50] <- 0
#' s <- 1
#' x <- theta + rnorm(200, 0, s)
#' pn.res <- ebnm_point_normal(x, s)
#' plot(pn.res)
#'
#' pe.res <- ebnm_point_exponential(x, s)
#' plot(pn.res, pe.res)
#'
#' # Customize plot:
#' library(ggplot2)
#' plot(pn.res, pe.res, remove_abline = TRUE) +
#'   theme_bw() +
#'   labs(x = "Simulated data")
#'
#' @export
#'
plot.ebnm <- function(x,
                      ...,
                      incl_pm = TRUE,
                      incl_cdf = FALSE,
                      subset = NULL,
                      remove_abline = FALSE) {
  # Hack to appease R CMD check ("no visible binding for global variable" note):
  obs <- pm <- label <- y <- NULL

  if (!inherits(x, "ebnm")) {
    stop("Input argument x must be an instance of class \"ebnm\".")
  }

  if (is.null(x[[data_ret_str()]])) {
    stop("Data not found in ebnm object. Results cannot be plotted.")
  }
  if (incl_pm && is.null(x[[df_ret_str()]][[pm_ret_str()]])) {
    warning("Posterior means not found in ebnm object and will not be plotted.")
    incl_pm <- FALSE
  }
  if (incl_cdf && is.null(x[[g_ret_str()]])) {
    warning("Fitted prior not found in ebnm object. CDF will not be plotted.")
    incl_cdf <- FALSE
  }

  ebnm_label <- deparse(substitute(x))
  llik <- as.numeric(x[[llik_ret_str()]])
  if (!is.null(llik)) {
    ebnm_label <- paste0(
      ebnm_label, " (llik: ", format(round(llik, 2), nsmall = 2), ")"
    )
  }

  if (incl_pm) {
    df1 <- data.frame(
      idx = 1:nrow(x[[data_ret_str()]]),
      name = rownames(x[[data_ret_str()]]),
      obs = x[[data_ret_str()]][[obs_ret_str()]],
      pm = x[[df_ret_str()]][[pm_ret_str()]],
      label = factor(ebnm_label)
    )
  }
  if (incl_cdf) {
    g_list <- list()
    g_list[[ebnm_label]] <- x[[g_ret_str()]]
  }

  args <- list(...)
  varnames <- sapply(substitute(list(...))[-1], deparse)
  if (length(args) > 0) {
    ebnm_idx <- which(sapply(args, inherits, "ebnm"))
    for (idx in ebnm_idx) {
      next_ebnm <- args[[idx]]
      if (is.null(next_ebnm[[data_ret_str()]])) {
        warning("An additional ebnm object was included as argument, but it does ",
                "not include a data field. Object will be ignored.")
      } else if (incl_pm && is.null(next_ebnm[[df_ret_str()]][[pm_ret_str()]])) {
        warning("An additional ebnm object was included as argument, but it does ",
                "not include posterior means. Object will be ignored.")
      } else if (incl_cdf && is.null(next_ebnm[[g_ret_str()]])) {
        warning("An additional ebnm object was included as argument, but it does ",
                "not include fitted prior. Object will be ignored.")
      } else if (!identical(x[[data_ret_str()]], next_ebnm[[data_ret_str()]])) {
        warning("An additional ebnm object was included as argument, but a different ",
                "dataset was used to fit the model. Object will be ignored.")
      } else {
        ebnm_label <- varnames[idx]
        llik <- as.numeric(next_ebnm[[llik_ret_str()]])
        if (!is.null(llik)) {
          ebnm_label <- paste0(
            ebnm_label, " (llik: ", format(round(llik, 2), nsmall = 2), ")"
          )
        }

        if (incl_pm) {
          next_df <- data.frame(
            idx = 1:nrow(next_ebnm[[data_ret_str()]]),
            name = rownames(next_ebnm[[data_ret_str()]]),
            obs = next_ebnm[[data_ret_str()]][[obs_ret_str()]],
            pm = next_ebnm[[df_ret_str()]][[pm_ret_str()]]
          )

          next_df$label <- ebnm_label

          levels(df1$label) <- c(levels(df1$label), ebnm_label)
          df1 <- rbind(df1, next_df)
        }
        if (incl_cdf) {
          g_list[[ebnm_label]] <- next_ebnm[[g_ret_str()]]
        }
      }
    }
    args <- args[-ebnm_idx]
    if (length(args) > 0) {
      warning("Additional arguments not of class ebnm were included. They will ",
              "be ignored.")
    }
  }

  if (incl_pm && is.numeric(subset)) {
    df1 <- df1[df1$idx %in% subset, ]
  } else if (incl_pm && is.character(subset)) {
    df1 <- df1[df1$name %in% subset, ]
  } else if (incl_pm && !is.null(subset)) {
    warning("Argument to subset must be NULL or a numeric or character vector. ",
            "It will be ignored.")
  }

  all_plots <- list()

  if (incl_pm) {
    if (length(unique(df1$label)) > 1) {
      p1 <- ggplot(df1, aes(x = obs, y = pm, color = label)) +
        geom_point() +
        scale_color_brewer(palette = "Set1") +
        labs(x = "Observations", y = "Posterior means",
             color = "Fitted EBNM models")
    } else {
      p1 <- ggplot(df1, aes(x = obs, y = pm)) +
        geom_point(color = "dodgerblue") +
        labs(x = "Observations", y = "Posterior means",
             title = paste("Log likelihood for model:", round(llik, 2)))
    }

    p1 <- p1 + theme_minimal()

    if (!remove_abline) {
      p1 <- p1 + geom_abline(slope = 1, linetype = "dashed")
    }

    all_plots[["pm"]] <- p1
  }

  if (incl_cdf) {
    if (incl_pm) {
      xgrid <- seq(min(df1$obs), max(df1$obs), length.out = 300)
    } else {
      xrange <- range(x[[data_ret_str()]][[obs_ret_str()]])
      xgrid <- seq(min(xrange), max(xrange), length.out = 300)
    }

    df2 <- data.frame(
      x = rep(xgrid, times = length(g_list)),
      y = as.vector(sapply(1:length(g_list), function(i) mixcdf(g_list[[i]], xgrid))),
      label = rep(factor(names(g_list), levels = names(g_list)), each = length(xgrid))
    )

    if (length(unique(df2$label)) > 1) {
      p2 <- ggplot(df2, aes(x = x, y = y, color = label)) +
        geom_line() +
        scale_color_brewer(palette = "Set1") +
        labs(x = expression(theta), y = "Cumulative prior probability",
             color = "Fitted EBNM models", title = "CDFs of fitted priors")
    } else {
      p2 <- ggplot(df2, aes(x = x, y = y)) +
        geom_line(color = "dodgerblue") +
        labs(x = expression(theta), y = "Cumulative prior probability",
             title = "CDF of fitted prior")
    }

    p2 <- p2 + theme_minimal()

    all_plots[["cdf"]] <- p2
  }

  if (length(all_plots) == 0) {
    cat("Nothing to plot.\n")
    return(NULL)
  } else if (length(all_plots) == 1) {
    return(all_plots[[1]])
  } else {
    oask <- devAskNewPage(TRUE)
    print(all_plots[["cdf"]])
    print(all_plots[["pm"]])
    devAskNewPage(oask)
    return(invisible())
  }
}

#' Summarize an ebnm object
#'
#' The \code{summary} method for class \code{\link{ebnm}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @return A \code{summary.ebnm} object.
#'
#' @method summary ebnm
#'
#' @export
#'
summary.ebnm <- function(object, ...) {
  retlist <- list()

  dat <- object[[data_ret_str()]]
  if (!is.null(dat)) {
    retlist$nobs <- length(dat[[obs_ret_str()]])
    retlist$heteroskedastic <- (diff(range(dat[[se_ret_str()]])) > 0)
  } else {
    retlist$nobs <- NA
    retlist$heteroskedastic <- NA
  }

  g <- object[[g_ret_str()]]

  if (!is.null(g)) {
    retlist$prior_family <- infer_prior_family(g)

    # Identify pointmass.
    if (inherits(g, "normalmix")) {
      pointmass_idx <- which(g$sd == 0)
    } else if (inherits(g, c("laplacemix", "gammamix"))) {
      pointmass_idx <- which(g$scale == 0)
    } else if (inherits(g, "unimix")) {
      pointmass_idx <- which(g$a == g$b)
    } else {
      pointmass_idx <- numeric(0)
    }
    if (length(pointmass_idx) == 1) {
      if (inherits(g, c("normalmix", "laplacemix"))) {
        retlist$pointmass_location <- g$mean[pointmass_idx]
      } else if (inherits(g, "gammamix")) {
        retlist$pointmass_location <- g$shift[pointmass_idx]
      } else if (inherits(g, "unimix")) {
        retlist$pointmass_location <- g$a[pointmass_idx]
      } else {
        retlist$pointmass_location <- NA
      }
      retlist$pointmass_weight <- g$pi[pointmass_idx]
    } else {
      retlist$pointmass_weight <- NA
    }
  } else {
    retlist$prior_family <- NA
    retlist$pointmass_location <- NA
    retlist$pointmass_weight <- NA
  }

  retlist$log_likelihood <- object[[llik_ret_str()]]
  if (is.null(retlist$log_likelihood)) {
    retlist$log_likelihood <- NA
  }

  retlist$posterior_summaries <- names(object[[df_ret_str()]])
  if (is.null(retlist$posterior_summaries)) {
    retlist$posterior_summaries <- NA
  }

  retlist$sampler_included <- !is.null(object[[samp_ret_str()]])

  retlist$call <- object[["call"]]

  class(retlist) <- c("summary.ebnm", "list")

  return(retlist)
}

#' Print a summary.ebnm object
#'
#' The \code{print} method for class \code{\link{summary.ebnm}}.
#'
#' @param x The \code{summary.ebnm} object.
#'
#' @param digits Number of significant digits to use.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @method print summary.ebnm
#'
#' @export
#'
print.summary.ebnm <- function(x, digits = 2, ...) {
  print_it(x, digits, summary = TRUE)
}

#' Print an ebnm object
#'
#' The \code{print} method for class \code{\link{ebnm}}.
#'
#' @param x The fitted \code{ebnm} object.
#'
#' @param digits Number of significant digits to use.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @method print ebnm
#'
#' @export
#'
print.ebnm <- function(x, digits = 2, ...) {
  print_it(summary(x), digits, summary = FALSE)
}

print_it <- function(x, digits, summary) {
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
  }

  if (!is.null(x$nobs)) {
    cat("EBNM model was fitted to", x$nobs, "observations with",
        ifelse(x$heteroskedastic, "_heteroskedastic_", "_homoskedastic_"),
        "standard errors.\n")
    cat("\n")
  }

  if (!is.null(x$prior_family) && x$prior_family != "npmle") {
    if (x$prior_family == "unknown") {
      cat("The fitted prior belongs to an _unrecognized_ prior family.\n")
    } else {
      cat("The fitted prior belongs to the",
          paste0("_", x$prior_family, "_"),
          "prior family.\n")
    }
    if (!is.null(x$pointmass_location)) {
      cat("It includes a point mass at",
          signif(x$pointmass_location, digits = digits),
          "with component weight equal to",
          paste0(signif(x$pointmass_weight, digits = digits), ".\n"))
    }
    cat("\n")
  }

  if (!is.null(x$log_likelihood)) {
    if (!is.null(attr(x$log_likelihood, "df"))) {
      cat(attr(x$log_likelihood, "df"),
          "degrees of freedom were used to estimate the model.\n")
    }
    cat("The log likelihood for the model is",
        paste0(round(as.numeric(x$log_likelihood), digits = digits), ".\n\n"))
  }

  if (summary) {
    if (length(x$posterior_summaries) > 0) {
      cat("Available posterior summaries:",
          paste0("_", paste(x$posterior_summaries, collapse = "_, _"), "_.\n"))
      if (pm_ret_str() %in% x$posterior_summaries) {
        cat("Use method fitted() to access available summaries.\n")
      }
    } else {
      cat("Posterior summaries are not available.\n")
    }
    cat("\n")

    if (x$sampler_included) {
      cat("A posterior sampler _is_ available and can be accessed using",
          "method simulate().\n")
    } else {
      cat("A posterior sampler is _not_ available.\n")
      if (x$prior_family != "unknown") {
        cat("One can be added via function ebnm_add_sampler().\n")
      }
    }
    cat("\n")
  }

  return(invisible(x))
}

#' Extract the log likelihood from a fitted EBNM model
#'
#' The \code{\link[stats]{logLik}} method for class \code{\link{ebnm}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @return An object of class \code{logLik}, which includes attributes \code{df},
#'   the degrees of freedom --- i.e., number of parameters in the model ---, and
#'   \code{nobs}, the number of observations in the data.
#'
#' @importFrom stats logLik
#' @method logLik ebnm
#'
#' @export
#'
logLik.ebnm <- function(object, ...) {
  return(object[[llik_ret_str()]])
}

#' Extract posterior estimates from a fitted EBNM model
#'
#' The \code{\link[stats]{fitted}} method for class \code{\link{ebnm}}.
#'   Returns a data frame that includes posterior means, standard deviations,
#'   and local false sign rates (when available).
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @importFrom stats fitted
#' @method fitted ebnm
#'
#' @export
#'
fitted.ebnm <- function(object, ...) {
  return(object[[df_ret_str()]])
}

#' Extract posterior means from a fitted EBNM model
#'
#' The \code{\link[stats]{coef}} method for class \code{\link{ebnm}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @importFrom stats coef
#' @method coef ebnm
#'
#' @export
#'
coef.ebnm <- function(object, ...) {
  posterior <- fitted(object)
  retval <- posterior[[pm_ret_str()]]
  names(retval) <- rownames(posterior)
  return(retval)
}

#' Extract posterior variances from a fitted EBNM model
#'
#' The \code{\link[stats]{vcov}} method for class \code{\link{ebnm}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @importFrom stats vcov
#' @method vcov ebnm
#'
#' @export
#'
vcov.ebnm <- function(object, ...) {
  posterior <- fitted(object)
  retval <- posterior[[psd_ret_str()]]^2
  names(retval) <- rownames(posterior)
  return(retval)
}

#' Calculate residuals for a fitted EBNM model
#'
#' The \code{\link[stats]{residuals}} method for class \code{\link{ebnm}}.
#'   Calculates "residuals" \eqn{x_i - \hat{\theta_i}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @importFrom stats residuals
#' @method residuals ebnm
#'
#' @export
#'
residuals.ebnm <- function(object, ...) {
  retval <- object[[data_ret_str()]][[obs_ret_str()]] -
    object[[df_ret_str()]][[pm_ret_str()]]
  names(retval) <- rownames(object[[data_ret_str()]])
  return(retval)
}

#' Use the estimated prior from a fitted EBNM model to solve the EBNM problem for
#'   new data
#'
#' The \code{\link[stats]{predict}} method for class \code{\link{ebnm}}.
#'
#' @inheritParams ebnm
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param newdata A vector of new observations. Missing observations
#'   (\code{NA}s) are not allowed.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @return A data frame that includes posterior means, posterior standard
#'   deviations, and local false sign rates for the observations in \code{newdata}.
#'
#' @importFrom stats predict
#' @method predict ebnm
#'
#' @export
#'
predict.ebnm <- function(object, newdata, s = 1, ...) {
  g_init <- object[[g_ret_str()]]
  ebnm_res <- ebnm(
    newdata,
    s,
    g_init = g_init,
    fix_g = TRUE,
    output = c(pm_arg_str(), psd_arg_str(), lfsr_arg_str())
  )
  return(fitted(ebnm_res))
}

#' Get the number of observations used to fit an EBNM model
#'
#' The \code{\link[stats]{nobs}} method for class \code{\link{ebnm}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @return The number of observations used to fit the \code{ebnm} object.
#'
#' @importFrom stats nobs
#' @method nobs ebnm
#'
#' @export
#'
nobs.ebnm <- function(object, ...) {
  retval <- length(object[[data_ret_str()]][[obs_ret_str()]])
  if (retval > 0) {
    return(retval)
  } else {
    return(NULL)
  }
}

#' Sample from the posterior of a fitted EBNM model
#'
#' The \code{\link[stats]{simulate}} method for class \code{\link{ebnm}}.
#'   Extracts the posterior sampler from an object of class \code{\link{ebnm}}
#'   and returns a specified number of samples.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param nsim The number of posterior samples to return per observation.
#'
#' @param seed Either \code{NULL} or an integer that will be used in a call to
#'   \code{set.seed} before simulating. If set, the value is saved as the
#'   \code{"seed"} attribute of the returned value. The default, \code{NULL},
#'   will not change the random generator state.
#'
#' @param ... Additional arguments to be passed to the posterior sampler
#'   function. Since \code{ebnm_horseshoe} returns an MCMC sampler, it takes
#'   parameter \code{burn}, the number of burn-in samples to discard.  At
#'   present, no other samplers take any additional parameters.
#'
#' @return A matrix of posterior samples, with rows corresponding to
#'   distinct samples and columns corresponding to observations.
#'
#' @importFrom stats simulate
#' @method simulate ebnm
#'
#' @export
#'
simulate.ebnm <- function(object, nsim = 1, seed = NULL, ...) {
  if (is.null(object[[samp_ret_str()]])) {
    stop("Object does not contain a posterior sampler. Note that samplers are ",
         "not returned by default. One can be added via function ",
         "ebnm_add_sampler().")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  retval <- object[[samp_ret_str()]](nsim, ...)
  if (!is.null(seed)) {
    attr(retval, "seed") <- seed
  }
  return(retval)
}

#' Obtain posterior quantiles using a fitted EBNM model
#'
#' The \code{\link[stats]{quantile}} method for class \code{\link{ebnm}}.
#'   Quantiles for posterior distributions \eqn{\theta_i \mid x_i, s_i, g} are
#'   estimated via sampling. By default, \code{\link{ebnm}} does not return a
#'   posterior sampler; one can be added to the \code{ebnm} object using
#'   function \code{\link{ebnm_add_sampler}}.
#'
#' @inheritParams stats::quantile
#'
#' @param x The fitted \code{ebnm} object.
#'
#' @param type An integer between 1 and 9 selecting one of the nine quantile
#'   algorithms detailed in \code{\link[stats]{quantile}} to be used.
#'
#' @param nsim The number of samples to use to estimate quantiles.
#'
#' @param ... Additional arguments to be passed to the posterior sampler
#'   function. Since \code{ebnm_horseshoe} returns an MCMC sampler, it takes
#'   parameter \code{burn}, the number of burn-in samples to discard.  At
#'   present, no other samplers take any additional parameters.
#'
#' @return A matrix with columns giving quantiles for each posterior
#'   \eqn{\theta_i \mid x_i, s_i, g}.
#'
#' @importFrom stats quantile
#' @method quantile ebnm
#'
#' @export
#'
quantile.ebnm <- function(x, probs = seq(0, 1, 0.25),
                          names = TRUE, type = 7, digits = 7, nsim = 1000, ...) {
  if (is.null(x[[samp_ret_str()]])) {
    stop("Quantiles are estimated by sampling from the posterior. Note that ",
         "samplers are not returned by default. One can be added via ",
         "function ebnm_add_sampler().")
  }
  samp <- simulate(x, nsim = nsim, ...)
  return(t(apply(samp, 2, quantile, probs = probs,
                 names = names, type = type, digits = digits)))
}

#' Obtain credible intervals using a fitted EBNM model
#'
#' The \code{\link[stats]{confint}} method for class \code{\link{ebnm}}.
#'   Estimates posterior "credible intervals" for each "true mean" \eqn{\theta_i}.
#'   We define the \eqn{(1 - \alpha)}\% credible interval for \eqn{\theta_i} as
#'   the narrowest continuous interval \eqn{[a_i, b_i]} such that
#'   \eqn{\theta_i \in [a_i, b_i]} with posterior probability at least
#'   \eqn{1 - \alpha}, where \eqn{\alpha \in (0,1)}. We estimate these credible
#'   intervals using Monte Carlo sampling. Note
#'   that by default, \code{\link{ebnm}} does not return a posterior
#'   sampler; one can be added to the \code{ebnm} object using function
#'   \code{\link{ebnm_add_sampler}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param parm A vector of numeric indices specifying which means \eqn{\theta_i}
#'   are to be given confidence intervals. If missing, all observations are
#'   considered.
#'
#' @param level The "confidence level" \eqn{1 - \alpha} desired.
#'
#' @param nsim The number of samples to use to estimate confidence intervals.
#'
#' @param ... Additional arguments to be passed to the posterior sampler
#'   function. Since \code{ebnm_horseshoe} returns an MCMC sampler, it takes
#'   parameter \code{burn}, the number of burn-in samples to discard.  At
#'   present, no other samplers take any additional parameters.
#'
#' @return A matrix with columns giving lower and upper confidence limits for
#'   each mean \eqn{\theta_i}. These will be labelled as "CI.lower" and
#'   "CI.upper".
#'
#' @importFrom stats confint
#' @method confint ebnm
#'
#' @export
#'
confint.ebnm <- function(object, parm, level = 0.95, nsim = 1000, ...) {
  if (is.null(object[[samp_ret_str()]])) {
    stop("Confidence intervals are obtained by sampling from the posterior. ",
         "Note that samplers are not returned by default. One can be added ",
         "via function ebnm_add_sampler().")
  }

  samp <- simulate(object, nsim = nsim, ...)
  if (!missing(parm)) {
    samp <- samp[, parm, drop = FALSE]
  }
  samp <- apply(samp, 2, sort)

  m <- round(nsim * (1 - level))
  y <- apply(samp, 2, function(x) x[seq(nsim - m + 1, nsim)] - x[seq(1, m)])
  i <- apply(y, 2, which.min)
  hpd <- t(sapply(1:ncol(samp),
                  function(j) c(samp[i[j], j], samp[nsim - m + i[j], j])))
  rownames(hpd) <- colnames(samp)
  colnames(hpd) <- c("CI.lower", "CI.upper")

  return(hpd)
}
