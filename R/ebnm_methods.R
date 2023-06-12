#' Plot an ebnm object
#'
#' Given one or more fitted \code{\link{ebnm}} object(s), produces a plot of
#'   posterior means vs. observations.
#'
#' An object of class \code{ggplot} is returned, so that the plot can be
#'   customized in the usual \code{\link[ggplot2]{ggplot2}} fashion.
#'
#' @param x The fitted \code{ebnm} object.
#'
#' @param remove_abline To better illustrate shrinkage effects, the plot
#'   includes the line \eqn{y = x} by default. If
#'   \code{remove_abline = TRUE}, then this line will not be drawn.
#'
#' @param ... Additional \code{ebnm} objects to be included on the same plot.
#'
#' @return A \code{ggplot} object.
#'
#' @method plot ebnm
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom methods is
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
plot.ebnm <- function(x, ..., remove_abline = FALSE) {
  if (!inherits(x, "ebnm")) {
    stop("Input argument x must be an instance of class \"ebnm\".")
  }

  if (is.null(x[[data_ret_str()]])) {
    stop("Data not found in ebnm object. Results cannot be plotted.")
  }

  if (is.null(x[[df_ret_str()]][[pm_ret_str()]])) {
    stop("Posterior means not found in ebnm object. Results cannot be plotted.")
  }

  ebnm_label <- deparse(substitute(x))
  llik <- as.numeric(x[[llik_ret_str()]])
  if (!is.null(llik)) {
    ebnm_label <- paste0(
      ebnm_label, " (llik: ", format(round(llik, 2), nsmall = 2), ")"
    )
  }

  df <- data.frame(
    obs = x[[data_ret_str()]][[obs_ret_str()]],
    pm = x[[df_ret_str()]][[pm_ret_str()]],
    label = factor(ebnm_label)
  )

  args <- list(...)
  varnames <- sapply(substitute(list(...))[-1], deparse)
  if (length(args) > 0) {
    ebnm_idx <- which(sapply(args, inherits, "ebnm"))
    for (idx in ebnm_idx) {
      next_ebnm <- args[[idx]]
      if (is.null(next_ebnm[[data_ret_str()]])) {
        warning("An additional ebnm object was included as argument, but it does ",
                "not include a data field. Object will be ignored.")
      } else if (is.null(next_ebnm[[df_ret_str()]][[pm_ret_str()]])) {
        warning("An additional ebnm object was included as argument, but it does ",
                "not include posterior means. Object will be ignored.")
      } else if (!identical(x[[data_ret_str()]], next_ebnm[[data_ret_str()]])) {
        warning("An additional ebnm object was included as argument, but a different ",
                "dataset was used to fit the model. Object will be ignored.")
      } else {
        next_df <- data.frame(
          obs = next_ebnm[[data_ret_str()]][[obs_ret_str()]],
          pm = next_ebnm[[df_ret_str()]][[pm_ret_str()]]
        )

        ebnm_label <- varnames[idx]
        llik <- as.numeric(next_ebnm[[llik_ret_str()]])
        if (!is.null(llik)) {
          ebnm_label <- paste0(
            ebnm_label, " (llik: ", format(round(llik, 2), nsmall = 2), ")"
          )
        }
        next_df$label <- ebnm_label

        levels(df$label) <- c(levels(df$label), ebnm_label)
        df <- rbind(df, next_df)
      }
    }
    args <- args[-ebnm_idx]
    if (length(args) > 0) {
      warning("Additional arguments not of class ebnm were included. They will ",
              "be ignored.")
    }
  }

  if (length(unique(df$label)) > 1) {
    plt <- ggplot(df, aes(x = obs, y = pm, color = label)) +
      geom_point() +
      labs(x = "Observations", y = "Posterior means",
           color = "Fitted EBNM models") +
      theme_minimal()
  } else {
    plt <- ggplot(df, aes(x = obs, y = pm)) +
      geom_point() +
      labs(x = "Observations", y = "Posterior means",
           title = paste("Log likelihood for model:", round(llik, 2))) +
      theme_minimal()
  }

  if (!remove_abline) {
    plt <- plt + geom_abline(slope = 1, linetype = "dashed")
  }

  return(plt)
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
    retlist$prior_class <- class(g)

    # Identify pointmass.
    if (class(g) == "normalmix") {
      pointmass_idx <- which(g$sd == 0)
    } else if (class(g) %in% c("laplacemix", "gammamix")) {
      pointmass_idx <- which(g$scale == 0)
    } else if (class(g) == "unimix") {
      pointmass_idx <- which(g$a == g$b)
    } else {
      pointmass_idx <- numeric(0)
    }
    if (length(pointmass_idx) == 1) {
      if (class(g) %in% c("normalmix", "laplacemix")) {
        retlist$pointmass_location <- g$mean[pointmass_idx]
      } else if (class(g) == "gammamix") {
        retlist$pointmass_location <- g$shift[pointmass_idx]
      } else if (class(g) == "unimix") {
        retlist$pointmass_location <- g$a[pointmass_idx]
      } else {
        retlist$pointmass_location <- NA
      }
      retlist$pointmass_weight <- g$pi[pointmass_idx]
    } else {
      retlist$pointmass_weight <- NA
    }
  } else {
    retlist$prior_class <- NA
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
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  if (!is.null(x$nobs)) {
    cat("EBNM model was fitted to", x$nobs, "observations with",
        ifelse(x$heteroskedastic, "_heteroskedastic_", "_homoskedastic_"),
        "standard errors.\n")
    cat("\n")
  }

  if (!is.null(x$prior_class)) {
    cat("The fitted prior is of class", paste0("_", x$prior_class, "_.\n"))
    if (!is.null(x$pointmass_location)) {
      cat("It includes a point mass at",
          signif(x$pointmass_location, digits = digits),
          "with component weight equal to",
          paste0(signif(x$pointmass_weight, digits = digits), ".\n"))
    }
    cat("\n")
  }

  if (!is.null(x$log_likelihood)) {
    cat(attr(x$log_likelihood, "df"),
        "degrees of freedom were used to estimate the model.\n")
    cat("The likelihood is",
        paste0(round(as.numeric(x$log_likelihood), digits = digits), ".\n\n"))
  }

  if (summary) {
    if (length(x$posterior_summaries) > 0) {
      cat("Available posterior summaries:",
          paste0("_", paste(x$posterior_summaries, collapse = "_, _"), "_.\n"))
      if (pm_ret_str() %in% x$posterior_summaries) {
        cat("Use method fitted() with 'only_means = FALSE' to access all",
            "available summaries.\n")
      }
    } else {
      cat("Posterior summaries are not available.\n")
    }
    cat("\n")

    if (x$sampler_included) {
      cat("A posterior sampler _is_ available and can be accessed using",
          "method samp().\n")
    } else {
      cat("A posterior sampler is _not_ available.\nOne can be obtained by",
          "re-running ebnm() with argument 'output = output_all()'.\n")
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
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param only_means Return only posterior means, or return a data frame that
#'   includes standard deviations and local false sign rates (when available)?
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @method fitted ebnm
#'
#' @export
#'
fitted.ebnm <- function(object, only_means = TRUE, ...) {
  if (only_means) {
    return(object[[df_ret_str()]][[pm_ret_str()]])
  } else {
    return(object[[df_ret_str()]])
  }
}

#' Use the estimated prior from a fitted EBNM model to solve the EBNM problem for
#'   new data
#'
#' The \code{\link[stats]{predict}} method for class \code{\link{ebnm}}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param newdata A list that includes fields \code{x} (the new observations)
#'   and \code{s} (corresponding standard errors). See \code{\link{ebnm}} for
#'   details about how to specify \code{x} and \code{s}.
#'
#' @param output A character vector indicating which values are to be returned.
#'   See \code{\link{ebnm}} for details.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @return An \code{ebnm} object. See \code{\link{ebnm}} for details.
#'
#' @method predict ebnm
#'
#' @export
#'
predict.ebnm <- function(object, newdata, output = output_default(), ...) {
  if (!("x" %in% names(newdata) && "s" %in% names(newdata))) {
    stop("Argument 'newdata' must be a list with fields 'x' (the new ",
         "observations) and 's' (corresponding standard errors).")
  }
  g_init <- object[[g_ret_str()]]
  return(ebnm(
    newdata$x,
    newdata$s,
    g_init = g_init,
    fix_g = TRUE,
    output = output
  ))
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
#' A convenience function that extracts the posterior sampler from an object
#'   of class \code{\link{ebnm}} and returns a specified number of samples.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param nsamp The number of posterior samples to return per observation.
#'
#' @param ... Additional arguments to be passed to the posterior sampler
#'   function. Since \code{ebnm_horseshoe} returns an MCMC sampler, it takes
#'   parameter \code{burn}, the number of burn-in samples to discard.  At
#'   present, no other samplers take any additional parameters.
#'
#' @return A matrix of posterior samples, with rows corresponding to
#'   distinct samples and columns corresponding to observations.
#'
#' @method samp ebnm
#'
#' @export
#'
samp.ebnm <- function(object, nsamp, ...) {
  if (is.null(object[[samp_ret_str()]])) {
    stop("Object does not contain a posterior sampler. Note that samplers are ",
         "not returned by default. To obtain one, include argument ",
         "'output = output_all()' in the call to ebnm.")
  }
  object[[samp_ret_str()]](nsamp, ...)
}

samp <- function(object, ...) {
  UseMethod("samp")
}

#' Obtain confidence intervals using a fitted EBNM model
#'
#' The \code{\link[stats]{confint}} method for class \code{\link{ebnm}}.
#'   Confidence intervals for means \eqn{\theta_i} are obtained by sampling from
#'   the posterior. By default, \code{\link{ebnm}} does not return a posterior
#'   sampler; one can be obtained by setting \code{output = output_all()} in the
#'   call to \code{ebnm}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param level The confidence level required.
#'
#' @param nsamp The number of samples to use to estimate confidence intervals.
#'
#' @param ... Additional arguments to be passed to the posterior sampler
#'   function. Since \code{ebnm_horseshoe} returns an MCMC sampler, it takes
#'   parameter \code{burn}, the number of burn-in samples to discard.  At
#'   present, no other samplers take any additional parameters.
#'
#' @return A matrix with columns giving lower and upper confidence limits for
#'   each mean \eqn{\theta_i}.
#'
#' @method confint ebnm
#'
#' @export
#'
confint.ebnm <- function(object, level = 0.95, nsamp = 1000, ...) {
  if (is.null(object[[samp_ret_str()]])) {
    stop("Confidence intervals are obtained by sampling from the posterior. ",
         "To obtain a posterior sampler, include argument 'output = output_all()' ",
         "in the call to ebnm.")
  }
  samp <- samp(object, nsamp = nsamp, ...)
  return(t(apply(samp, 2, quantile, probs = 0.5 + c(-1, 1) * level / 2)))
}
