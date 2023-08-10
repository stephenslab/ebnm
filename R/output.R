#' @describeIn ebnm Lists the default return values.
#'
#' @export
#'
ebnm_output_default <- function() {
  return(c(data_arg_str(), pm_arg_str(), psd_arg_str(),
           g_arg_str(), llik_arg_str()))
}

#' @describeIn ebnm Lists all valid return values.
#'
#' @export
#'
ebnm_output_all <- function() {
  return(c(data_arg_str(), pm_arg_str(), psd_arg_str(), pm2_arg_str(),
           lfsr_arg_str(), g_arg_str(), llik_arg_str(), samp_arg_str()))
}

# Return value names as used in arguments to parameter 'output'.
data_arg_str <- function() "data"
pm_arg_str   <- function() "posterior_mean"
psd_arg_str  <- function() "posterior_sd"
pm2_arg_str  <- function() "posterior_second_moment"
lfsr_arg_str <- function() "lfsr"
g_arg_str    <- function() "fitted_g"
llik_arg_str <- function() "log_likelihood"
samp_arg_str <- function() "posterior_sampler"

# Return value names as used in the returned ebnm object.
data_ret_str <- function() "data"
obs_ret_str  <- function() "x"
se_ret_str   <- function() "s"
df_ret_str   <- function() "posterior"
pm_ret_str   <- function() "mean"
psd_ret_str  <- function() "sd"
pm2_ret_str  <- function() "second_moment"
lfsr_ret_str <- function() "lfsr"
g_ret_str    <- function() "fitted_g"
llik_ret_str <- function() "log_likelihood"
samp_ret_str <- function() "posterior_sampler"
grp_ret_str  <- function() "group"

# Postprocessing of the returned object is done here.
as_ebnm <- function(retlist, call) {
  retlist$call <- call
  class(retlist) <- c("ebnm", "list")
  return(retlist)
}

ash_output <- function(output) {
  ash_arg_str <- c("data", "PosteriorMean", "PosteriorSD", "PosteriorSD",
                   "lfsr", "fitted_g", "loglik", "post_sampler")
  which_args  <- pmatch(output, ebnm_output_all())
  return(ash_arg_str[which_args])
}

posterior_in_output <- function(output) {
  post_args <- c(pm_arg_str(), psd_arg_str(), pm2_arg_str(), lfsr_arg_str())
  return(any(post_args %in% output))
}

result_in_output <- function(output) {
  res_args <- c(pm_arg_str(), psd_arg_str(), pm2_arg_str())
  return(any(res_args %in% output))
}

lfsr_in_output <- function(output) {
  return(lfsr_arg_str() %in% output)
}

data_in_output <- function(output) {
  return(data_arg_str() %in% output)
}

add_data_to_retlist <- function(retlist, x, s) {
  df <- list()
  df[[obs_ret_str()]] <- x
  df[[se_ret_str()]] <- rep(s, length.out = length(x))
  df <- data.frame(df)

  retlist[[data_ret_str()]] <- df
  return(retlist)
}

add_posterior_to_retlist <- function(retlist, posterior, output, x) {
  if (!posterior_in_output(output)) {
    return(retlist)
  }

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
  df <- data.frame(df, row.names = names(x))

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

add_llik_to_retlist <- function(retlist, llik, x, df) {
  attr(llik, "nobs") <- length(x)
  attr(llik, "df") <- df
  class(llik) <- "logLik"
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
