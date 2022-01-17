#' Plot an ebnm object
#'
#' Given a fitted ebnm object, produces a plot of posterior means vs.
#'   observations.
#'
#' An object of class \code{ggplot} is returned, so that the plot can be
#'   customized in the usual \code{\link[ggplot2]{ggplot2}} fashion.
#'
#' @param x The fitted ebnm object.
#'
#' @param remove_abline To better illustrate shrinkage effects, the plot
#'   will include the line \eqn{y = x} by default. If
#'   \code{remove_abline = TRUE}, then this line will not be drawn.
#'
#' @param ... Additional parameters to be passed to \code{ggplot2} function
#'   \code{\link[ggplot2]{geom_point}}.
#'
#' @method plot ebnm
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs theme_minimal
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

  if (is.null(x$data)) {
    stop("Data not found in ebnm object. Results cannot be plotted.")
  }

  if (is.null(x$posterior$mean)) {
    stop("Posterior means not found in ebnm object. Results cannot be plotted.")
  }

  df <- data.frame(
    x = x$data$x,
    pm = x$posterior$mean
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
