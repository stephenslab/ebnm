#' Add sampler to an ebnm_object
#'
#' Adds a posterior sampler to a fitted \code{\link{ebnm}} object.
#'
#' @param ebnm_res The fitted \code{ebnm} object.
#'
#' @return The \code{ebnm} object with an additional field
#'   \code{posterior_sampler}.
#'
#' @export
#'
ebnm_add_sampler <- function(ebnm_res) {
  if (!inherits(ebnm_res, "ebnm")) {
    stop("Input argument must be an instance of class \"ebnm\".")
  }
  if (is.null(ebnm_res[[data_ret_str()]])) {
    stop("Data not found in ebnm object. Sampler cannot be added.")
  }
  if (is.null(ebnm_res[[g_ret_str()]])) {
    warning("Fitted prior not found in ebnm object. Sampler cannot be added.")
    incl_cdf <- FALSE
  }

  x <- ebnm_res[[data_ret_str()]][[obs_ret_str()]]
  names(x) <- rownames(ebnm_res[[data_ret_str()]])
  sampler <- ebnm(
    x = x,
    s = ebnm_res[[data_ret_str()]][[se_ret_str()]],
    g_init = ebnm_res[[g_ret_str()]],
    fix_g = TRUE,
    output = samp_arg_str()
  )
  ebnm_res[[samp_ret_str()]] <- sampler[[samp_ret_str()]]
  return(ebnm_res)
}
