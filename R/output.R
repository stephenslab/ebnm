#' @describeIn ebnm Defines the default return values.
#'
#' @export
#'
output_default <- function() {
  return(c(pm_arg_str(), psd_arg_str(), g_arg_str(), llik_arg_str()))
}

#' @describeIn ebnm Lists all valid return values.
#'
#' @export
#'
output_all <- function() {
  return(c(pm_arg_str(), psd_arg_str(), pm2_arg_str(), lfsr_arg_str(),
           g_arg_str(), llik_arg_str(), samp_arg_str()))
}

# Return value names as used in arguments to parameter 'output'.
pm_arg_str   <- function() "posterior_mean"
psd_arg_str  <- function() "posterior_sd"
pm2_arg_str  <- function() "posterior_second_moment"
lfsr_arg_str <- function() "lfsr"
g_arg_str    <- function() "fitted_g"
llik_arg_str <- function() "log_likelihood"
samp_arg_str <- function() "posterior_sampler"

# Return value names as used in the returned ebnm object.
df_ret_str   <- function() "posterior"
pm_ret_str   <- function() "mean"
psd_ret_str  <- function() "sd"
pm2_ret_str  <- function() "second_moment"
lfsr_ret_str <- function() "lfsr"
g_ret_str    <- function() "fitted_g"
llik_ret_str <- function() "log_likelihood"
samp_ret_str <- function() "posterior_sampler"

# Postprocessing of the returned object is done here.
as_ebnm <- function(retlist) {
  class(retlist) <- c("ebnm", "list")
  return(retlist)
}

ash_output <- function(output) {
  ash_arg_str <- c("PosteriorMean", "PosteriorSD", "PosteriorSD", "lfsr",
                   "fitted_g", "loglik", "post_sampler")
  which_args  <- pmatch(output, output_all())
  return(ash_arg_str[which_args])
}

posterior_in_output <- function(output) {
  post_args <- c(pm_arg_str(), psd_arg_str(), pm2_arg_str(), lfsr_arg_str())
  return(any(post_args %in% output))
}

result_in_output <- function(output) {
  res_args <-c(pm_arg_str(), psd_arg_str(), pm2_arg_str())
  return(any(res_args %in% output))
}

lfsr_in_output <- function(output) {
  return(lfsr_arg_str() %in% output)
}

# 'posterior' should be a list with fields 'mean', 'sd', 'mean2', and 'lfsr'.
add_posterior_to_retlist <- function(retlist, posterior, output) {
  df <- list()
  if (pm_arg_str() %in% output) {
    df[[pm_ret_str()]] <- posterior$mean
  }
  if (psd_arg_str() %in% output) {
    df[[psd_ret_str()]] <- posterior$sd
  }
  if (pm2_arg_str() %in% output) {
    df[[pm2_ret_str()]] <- posterior$mean2
  }
  if (lfsr_arg_str() %in% output) {
    df[[lfsr_ret_str()]] <- posterior$lfsr
  }
  df <- data.frame(df)

  retlist[[df_ret_str()]] <- df
  return(retlist)
}

g_in_output <- function(output) {
  return(g_arg_str() %in% output)
}

add_g_to_retlist <- function(retlist, g) {
  retlist[[g_ret_str()]] <- g
  return(retlist)
}

llik_in_output <- function(output) {
  return(llik_arg_str() %in% output)
}

add_llik_to_retlist <- function(retlist, llik) {
  retlist[[llik_ret_str()]] <- llik
  return(retlist)
}

sampler_in_output <- function(output) {
  return(samp_arg_str() %in% output)
}

add_sampler_to_retlist <- function(retlist, sampler) {
  retlist[[samp_ret_str()]] <- sampler
  return(retlist)
}
