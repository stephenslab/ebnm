#' Plot an ebnm object
#'
#' Given a fitted \code{ebnm} object, produces a plot of posterior means vs.
#'   observations.
#'
#' An object of class \code{ggplot} is returned, so that the plot can be
#'   customized in the usual \code{\link[ggplot2]{ggplot2}} fashion.
#'
#' @param x The fitted \code{ebnm} object.
#'
#' @param remove_abline To better illustrate shrinkage effects, the plot
#'   will include the line \eqn{y = x} by default. If
#'   \code{remove_abline = TRUE}, then this line will not be drawn.
#'
#' @param ... Additional parameters to be passed to \code{ggplot2} function
#'   \code{\link[ggplot2]{geom_point}}.
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
#' s <- 1
#' x <- theta + rnorm(200, 0, s)
#' ebnm.res <- ebnm(x, s)
#' plot(ebnm.res)
#'
#' # Customize plot:
#' library(ggplot2)
#' plot(ebnm.res, color = "blue", remove_abline = TRUE) +
#'   theme_bw() +
#'   labs(x = "Simulated data")
#'
#' @export
#'
plot.ebnm <- function(x, remove_abline = FALSE, ...) {
  if (!is(x,"ebnm")) {
    stop("Input argument x must be an instance of class \"ebnm\".")
  }

  if (is.null(x[[data_ret_str()]])) {
    stop("Data not found in ebnm object. Results cannot be plotted.")
  }

  if (is.null(x[[df_ret_str()]][[pm_ret_str()]])) {
    stop("Posterior means not found in ebnm object. Results cannot be plotted.")
  }

  df <- data.frame(
    x = x[[data_ret_str()]][[obs_ret_str()]],
    pm = x[[df_ret_str()]][[pm_ret_str()]]
  )

  plt <- ggplot(df, aes(x = x, y = pm)) +
    geom_point(...) +
    labs(x = "Observations", y = "Posterior means") +
    theme_minimal()

  if (!remove_abline) {
    plt <- plt + geom_abline(slope = 1, linetype = "dashed")
  }

  return(plt)
}

#' Summarize an ebnm object
#'
#' The \code{summary} method for class \code{ebnm}.
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
#' The \code{print} method for class \code{summary.ebnm}.
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
#' The \code{print} method for class \code{ebnm}.
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
          paste0(paste(x$posterior_summaries, collapse = ", "), ".\n"))
    } else {
      cat("Posterior summaries are not available.\n")
    }
    if (x$sampler_included) {
      cat("A posterior sampler _is_ available.\n")
    } else {
      cat("A posterior sampler is _not_ available.\n")
    }
  }

  return(invisible(x))
}

#' Extract the log likelihood from a fitted EBNM model
#'
#' The \code{\link[stats]{logLik}} method for class \code{ebnm}.
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

#' Extract the posterior means from a fitted EBNM model
#'
#' The \code{\link[stats]{fitted}} method for class \code{ebnm}.
#'
#' @param object The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @return The posterior means \eqn{\hat{\theta}} in a fitted \code{ebnm} model.
#'
#' @method fitted ebnm
#'
#' @export
#'
fitted.ebnm <- function(object, ...) {
  return(object[[df_ret_str()]][[pm_ret_str()]])
}

#' Use the estimated prior from a fitted EBNM model to solve the EBNM problem for
#'   new data
#'
#' The \code{\link[stats]{predict}} method for class \code{ebnm}.
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
#' The \code{\link[stats]{nobs}} method for class \code{ebnm}.
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
#' A convenience function that extracts the \code{posterior_sampler} from an
#'   object of class \code{ebnm} and returns a specified number of samples.
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
#' @return The number of observations used to fit the \code{ebnm} object.
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
