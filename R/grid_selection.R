#' @importFrom ashr normalmix
#'
init_g_for_npmle <- function(x,
                             s,
                             scale = "estimate",
                             min_K = 3,
                             max_K = 300,
                             KLdiv_target = 1 / length(x),
                             force_pointmass = FALSE) {
  # Default is to use point masses:
  normalmix_sd <- 0

  if (identical(scale, "estimate")) {
    if (force_pointmass || (diff(range(x))) / max_K < 3 * min(s)) {
      # Use point-mass mixture.
      scale <- ebnm_scale_npmle(
        x, s, min_K, max_K, KLdiv_target, pointmass = TRUE
      )
    } else {
      # Use Gaussian mixture.
      scale <- ebnm_scale_npmle(
        x, s, min_K, max_K, KLdiv_target, pointmass = FALSE
      )
      normalmix_sd <- scale / 2
    }
  }

  grid <- seq(min(x), max(x), by = scale)
  ncomp <- length(grid)
  g_init <- normalmix(pi = rep(1 / ncomp, ncomp),
                      mean = grid,
                      sd = normalmix_sd)

  return(g_init)
}

#' Set scale parameter for NPMLE and deconvolveR prior family
#'
#' The default method for setting the \code{scale} parameter for functions
#'   \code{\link{ebnm_npmle}} and \code{\link{ebnm_deconvolver}}.
#'
#' @inheritParams ebnm_npmle
#'
#' @param min_K The minimum number of components \eqn{K} to include in the
#'   mixture of point masses used to approximate the nonparametric family of
#'   all distributions.
#'
#' @param max_K The maximum number of components \eqn{K} to include in the
#'   approximating mixture of point masses.
#'
#' @param KLdiv_target The desired bound \eqn{\kappa} on the KL-divergence from
#'   the solution obtained using the approximating mixture to the exact
#'   solution. More precisely, the scale parameter is set such that given
#'   the exact MLE
#'   \deqn{\hat{g} := \mathrm{argmax}_{g \in G} L(g),}
#'   where \eqn{G} is the full nonparametric family, and given the MLE for the
#'   approximating family \eqn{\tilde{G}}
#'   \deqn{\tilde{g} := \mathrm{argmax}_{g \in \tilde{G}} L(g),}
#'   we have that
#'   \deqn{\mathrm{KL}(\hat{g} \ast N(0, s^2) \mid \tilde{g} \ast N(0, s^2)) \le \kappa,}
#'   where \eqn{\ast \ N(0, s^2)} denotes convolution with the normal error
#'   distribution (the derivation of the bound assumes homoskedastic
#'   observations). For details, see \strong{References} below.
#'
#' @param pointmass When the range of the data is so large that \code{max_K}
#'   point masses cannot provide a good approximation to the family of all
#'   distributions, then \code{ebnm} will instead use a mixture of normal
#'   distributions, with the standard deviation of each component equal to
#'   \code{scale}\eqn{/ 2}. Setting \code{pointmass = FALSE} gives the
#'   default \code{scale} for this mixture of normal distributions.
#'
#'   To be exact, \code{ebnm} uses a mixture of normal distributions rather than
#'   a mixture of point masses when
#'   \deqn{\frac{\max(x) - \min(x)}{\min(s)} > 3 \ \mathrm{max}_K;} for a
#'   rationale, see \strong{References} below. Note however that \code{ebnm}
#'   only uses a mixture of normal distributions when \code{scale = "estimate"};
#'   if parameter \code{scale} is set manually, then a mixture of point masses
#'   will be used in all cases. To use a mixture of normal distributions with
#'   the scale set manually, an object created by the constructor function
#'   \code{\link[ashr]{normalmix}} must be provided as argument to parameter
#'   \code{g_init} in function \code{\link{ebnm_npmle}} or function
#'   \code{\link{ebnm_deconvolver}}.
#'
#' @references
#' Jason Willwerscheid (2021).
#'   \emph{Empirical Bayes Matrix Factorization: Methods and Applications}.
#'   University of Chicago, PhD dissertation.
#'
#' @export
#'
ebnm_scale_npmle <- function(x,
                             s,
                             min_K = 3,
                             max_K = 300,
                             KLdiv_target = 1 / length(x),
                             pointmass = TRUE) {
  if (pointmass) {
    scale <- min(s) * (64 * KLdiv_target)^(1 / 4)
  } else {
    scale <- 2 * min(s) * sqrt(exp(2 * KLdiv_target) - 1)
  }

  xrange <- diff(range(x))
  ncomp <- min(max(min_K - 1, ceiling(xrange / scale)), max_K - 1) + 1

  return(xrange / (ncomp - 1))
}

#' Set scale parameter for scale mixtures of normals
#'
#' The default method for setting the \code{scale} parameter for function
#'   \code{\link{ebnm_normal_scale_mixture}}.
#'
#' @inherit ebnm_scale_npmle
#'
#' @param mode A scalar specifying the mode of the prior \eqn{g}.
#'
#' @param min_K The minimum number of components \eqn{K} to include in the
#'   finite mixture of normal distributions used to approximate the
#'   nonparametric family of scale mixtures of normals.
#'
#' @param max_K The maximum number of components \eqn{K} to include in the
#'   approximating mixture of normal distributions.
#'
#' @importFrom stats approx
#'
#' @export
#'
ebnm_scale_normalmix <- function(x,
                                 s,
                                 mode = 0,
                                 min_K = 3,
                                 max_K = 300,
                                 KLdiv_target = 1 / length(x)) {
  max_x2 <- max((x - mode)^2)
  min_s2 <- min(s)^2

  max_mult <- (max(max_x2 / min_s2, 1) + 1)^(1 / (min_K - 1)) + 1e-8
  min_mult <- max_mult^((min_K - 1) / (max_K - 1)) - 1e-8

  grid_mult <- approx(
    y = smngrid$m,
    x = log(smngrid$KL),
    xout = log(KLdiv_target),
    rule = 2
  )[["y"]]

  grid_mult <- min(max(grid_mult, min_mult), max_mult)

  K <- ceiling(log(max(max_x2 / min_s2, 1) + 1, base = grid_mult)) + 1
  scale <- sqrt(min_s2 * (grid_mult^(0:(K - 1)) - 1))

  return(scale)
}

#' Set scale parameter for nonparametric unimodal prior families
#'
#' The default method for setting the \code{scale} parameter for functions
#'   \code{\link{ebnm_unimodal}}, \code{\link{ebnm_unimodal_symmetric}},
#'   \code{\link{ebnm_unimodal_nonnegative}}, and
#'   \code{\link{ebnm_unimodal_nonpositive}}.
#'
#' @inherit ebnm_scale_normalmix
#'
#' @param min_K The minimum number of components \eqn{K} to include in the
#'   finite mixture of uniform distributions used to approximate the
#'   nonparametric family of unimodal distributions.
#'
#' @param max_K The maximum number of components \eqn{K} to include in the
#'   approximating mixture of uniform distributions.
#'
#' @importFrom stats approx
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @importFrom dplyr filter arrange pull
#'
#' @export
#'
ebnm_scale_unimix <- function(x,
                              s,
                              mode = 0,
                              min_K = 3,
                              max_K = 300,
                              KLdiv_target = 1 / length(x)) {
  max_x <- max(abs(x - mode))
  min_s <- min(s)

  max_K <- min(max_K, max(symmunigrid$idx))

  # Trim grids to have desired number of components:
  grids <- symmunigrid %>%
    filter(.data[["idx"]] <= max_K) %>%
    arrange(.data[["idx"]])

  # Remove grids that are too fine to span the range of the data:
  KLs <- grids %>%
    filter(.data[["idx"]] == max_K) %>%
    filter(.data[["loc"]] > max_x / min_s) %>%
    pull(.data[["KL"]])

  # Remove grids that are too coarse to give the required number of components:
  KLs <- grids %>%
    filter(.data[["KL"]] %in% KLs,.data[["idx"]] == min_K) %>%
    filter(.data[["loc"]] < max_x / min_s) %>%
    pull(.data[["KL"]])

  if (length(KLs) == 0) {
    # Bail and use the scale mixture of normals grid:
    scale <- ebnm_scale_normalmix(
      x = (x - mode) / min_s,
      s = 1,
      mode = 0,
      min_K = max_K,
      max_K = max_K,
      KLdiv_target = KLdiv_target
    )
  } else if (KLdiv_target <= min(KLs)) {
    scale <- grids %>%
             filter(.data[["KL"]] == min(KLs)) %>%
             pull(.data[["loc"]])
  } else if (KLdiv_target >= max(KLs)) {
    scale <- grids %>%
             filter(.data[["KL"]] == max(KLs)) %>%
             pull(.data[["loc"]])
  } else {
    # Select the two grids that are nearest the target KL:
    lower_KL <- max(KLs[KLs <= KLdiv_target])
    upper_KL <- min(KLs[KLs > KLdiv_target])
    lower_grid <- grids %>%
                  filter(.data[["KL"]] == lower_KL) %>%
                  pull(.data[["loc"]])
    upper_grid <- grids %>%
                  filter(.data[["KL"]] == upper_KL) %>%
                  pull(.data[["loc"]])

    # Interpolate:
    prop <- (log(KLdiv_target) - log(lower_KL)) / (log(upper_KL) - log(lower_KL))
    scale <- lower_grid + prop * (upper_grid - lower_grid)
  }

  # Trim grid and re-scale:
  scale <- min_s * scale
  K <- min(sum(scale < max_x) + 1, max_K)
  scale <- scale[1:K]

  return(scale)
}

get_ashr_grid <- function(x,
                          s,
                          mode = 0,
                          grid_mult = sqrt(2)) {
  # Adapted from ashr:::autoselect.mixsd.
  sigmamin <- min(s[s > 0]) / 10
  sigmamax <- max(8 * sigmamin, 2 * sqrt(max((x - mode)^2 - s^2, 0)))
  npoint <- ceiling(log2(sigmamax / sigmamin) / log2(grid_mult))
  scale <- grid_mult^((-npoint):0) * sigmamax
  scale <- c(0, scale)

  return(scale)
}
