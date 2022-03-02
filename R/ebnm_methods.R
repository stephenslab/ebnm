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
#' @importFrom ggplot2 ggplot aes_string geom_point geom_abline labs
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
    "x" = x[[data_ret_str()]][[obs_ret_str()]],
    "pm" = x[[df_ret_str()]][[pm_ret_str()]]
  )

  plt <- ggplot(df, aes_string(x = "x", y = "pm")) +
    geom_point(...) +
    labs(x = "Observations", y = "Posterior means") +
    theme_minimal()

  if (!remove_abline) {
    plt <- plt + geom_abline(slope = 1, linetype = "dashed")
  }

  return(plt)
}

#' Print an ebnm object
#'
#' The \code{print} method for class \code{ebnm}.
#'
#' @param x The fitted \code{ebnm} object.
#'
#' @param ... Not used. Included for consistency as an S3 method.
#'
#' @method print ebnm
#'
#' @export
#'
print.ebnm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")

  skipline <- FALSE
  if (!is.null(x[[data_ret_str()]])) {
    cat("Sample size:", length(x[[data_ret_str()]][[1]]), "\n")
    skipline <- TRUE
  }
  if (!is.null(x[[g_ret_str()]])) {
    cat("Prior class:", class(x[[g_ret_str()]]), "\n")
    skipline <- TRUE
  }
  if(skipline) {
    cat("\n")
  }

  skipline <- FALSE
  if (!is.null(x[[llik_ret_str()]])) {
    cat("Log likelihood of fitted model:", x[[llik_ret_str()]], "\n")
    skipline <- TRUE
  }
  if (!is.null(x[[df_ret_str()]])) {
    cat("Posterior summaries included: ")
    cat(paste(names(x[[df_ret_str()]]), collapse = ", "), "\n")
    skipline <- TRUE
  }
  if (!is.null(x[[samp_ret_str()]])) {
    cat("Posterior sampler included. \n")
    skipline <- TRUE
  }
  if(skipline) {
    cat("\n")
  }

  return(invisible(x))
}
