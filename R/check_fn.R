#' Check a custom ebnm function
#'
#' Checks inputs and outputs of an ebnm-style function. Designed to
#'   troubleshoot custom functions, especially those that are intended for use
#'   in the \code{flashier} package (e.g., as argument to the \code{ebnm_fn}
#'   parameter in function \code{\link[flashier]{flash}}).
#'
#' @param fn The function to be checked.
#'
#' @param x A test set (vector) of observations.
#'
#' @param s A test set of standard errors. Typically, ebnm-style functions
#'   should be able to accept a vector of standard errors or a scalar if all
#'   standard errors are identical. This is not always the case; for example,
#'   the horseshoe prior family requires homoskedastic standard errors.
#'
#' @return Prints a success message and silently returns 1 if all checks pass.
#'   Otherwise the function errors out.
#'
#' @examples
#' ebnm_check_fn(ebnm_normal, x = rnorm(10, sd = 2), s = 1)
#'
#' @export
#'
ebnm_check_fn <- function(fn, x, s) {
  required_args <- c("x", "s", "g_init", "fix_g", "output")
  if (!all(required_args %in% names(as.list(args(fn))))) {
    stop("Function must accept the following arguments: ",
         toString(required_args))
  }

  ebnm_res <- fn(x = x, s = s, output = ebnm_output_all())
  if (!inherits(ebnm_res, "ebnm") || !inherits(ebnm_res, "list")) {
    stop("Function should return a named list of class 'ebnm'.")
  }

  required_ret_fields <- c(df_ret_str(), g_ret_str(), llik_ret_str())
  if (!all(required_ret_fields %in% names(ebnm_res))) {
    stop("Return object must include the following fields: ",
         toString(required_ret_fields))
  }

  llik <- ebnm_res[[llik_ret_str()]]
  if (!(is.numeric(llik) && length(llik) == 1)) {
    stop("Returned likelihood should be numeric (if likelihood-based ",
         "methods are not used, return 0 rather than NA).")
  }

  df <- ebnm_res[[df_ret_str()]]
  if (!inherits(df, "data.frame")) {
    stop("Posterior summaries should be returned as a data frame.")
  }

  required_post_cols <- c(pm_ret_str(), pm2_ret_str())
  if (!all(required_post_cols %in% names(df))) {
    stop("Posterior summaries should include the following columns: ",
         toString(required_post_cols))
  }

  ebnm_res2 <- fn(x = x + 1, s = s,
                  g_init = ebnm_res[[g_ret_str()]],
                  fix_g = TRUE,
                  output = ebnm_output_all())
  if (!identical(ebnm_res[[g_ret_str()]], ebnm_res2[[g_ret_str()]])) {
    stop("Argument fix_g appears not to be working correctly.")
  }

  if (samp_ret_str() %in% names(ebnm_res)) {
    nsamp <- ifelse(length(x) == 2, 3, 2)
    sampler <- ebnm_res[[samp_ret_str()]]
    samp <- sampler(nsamp)
    if (!is.matrix(samp) || nrow(samp) != nsamp) {
      stop("Sampler appears not to be working correctly; it should return a ",
           "matrix where each row corresponds to a different posterior sample ",
           "and each column corresponds to a different observation in the data ",
           "(i.e., a different 'true mean' theta_i.")
    }
  }

  message("Function has passed all checks.")
  return(invisible(1))
}
