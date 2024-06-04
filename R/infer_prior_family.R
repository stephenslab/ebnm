infer_prior_family <- function(g) {
  prior_family <- "unknown"
  if (inherits(g, "normalmix")) {
    if (length(g$sd) == 1) {
      if (g$sd == Inf) {
        prior_family <- "flat"
      } else if (g$sd == 0) {
        prior_family <- "point_mass"
      } else {
        prior_family <- "normal"
      }
    } else if (length(g$sd) == 2
               && g$sd[1] == 0
               && g$mean[1] == g$mean[2]) {
      prior_family <- "point_normal"
    } else if (all(g$mean == g$mean[1])) {
      prior_family <- "normal_scale_mixture"
    } else {
      prior_family <- "npmle"
    }
  } else if (inherits(g, "unimix")) {
    if (all((g$a + g$b) / 2 == (g$a[1] + g$b[1]) / 2)) {
      prior_family <- "unimodal_symmetric"
    } else if (all(g$a == g$a[1])) {
      prior_family <- "unimodal_nonnegative"
    } else if (all(g$b == g$b[1])) {
      prior_family <- "unimodal_nonpositive"
    } else if (all(g$a == g$a[1] | g$b == g$a[1])
               || all(g$a == g$b[1] | g$b == g$b[1])) {
      prior_family <- "unimodal"
    } else {
      prior_family <- "npmle"
    }
  } else if (inherits(g, "horseshoe")) {
    prior_family <- "horseshoe"
  } else if (inherits(g, "laplacemix")) {
    if (length(g$scale) == 2
        && g$scale[1] == 0
        && g$mean[1] == g$mean[2]) {
      prior_family <- "point_laplace"
    }
  } else if (inherits(g, "gammamix")) {
    if (length(g$scale) == 2
        && g$scale[1] == 0
        && g$shape[1] == 1
        && g$shape[2] == 1
        && g$shift[1] == g$shift[2])
      prior_family <- "point_exponential"
  } else if (inherits(g, "tnormalmix")) {
    if (length(g$mean) == 2
        && g$mean[1] == 0
        && g$sd[1] == 0
        && g$a[2] == 0
        && g$b[2] == Inf)
      prior_family <- "generalized_binary"
  }
  return(prior_family)
}
