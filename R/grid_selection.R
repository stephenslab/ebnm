#' @importFrom ashr normalmix
#'
init_g_for_npmle <- function(x,
                             s,
                             scale = "estimate",
                             min_K = 3,
                             max_K = 300,
                             KLdiv_target = 1 / length(x),
                             force_pointmass = FALSE) {
  xrange <- diff(range(x))

  # See ebnm paper for rationale.
  if (!identical(scale, "estimate")) {
    grid <- seq(min(x), max(x), by = scale)
    ncomp <- length(grid)

    g_init <- ashr::normalmix(pi = rep(1 / ncomp, ncomp),
                              mean = grid,
                              sd = 0)
  } else if (force_pointmass || xrange / max_K < 3 * min(s)) {
    # Use point-mass mixture.
    d <- min(s) * (64 * KLdiv_target)^(1 / 4)
    ncomp <- min(max(min_K - 1, ceiling(xrange / d)), max_K - 1) + 1
    grid <- seq(min(x), max(x), length.out = ncomp)

    g_init <- ashr::normalmix(pi = rep(1 / ncomp, ncomp),
                              mean = grid,
                              sd = 0)
  } else {
    # Use Gaussian mixture.
    d <- 2 * min(s) * sqrt(exp(2 * KLdiv_target) - 1)
    ncomp <- min(max(min_K - 1, ceiling(xrange / d)), max_K - 1) + 1
    grid <- seq(min(x), max(x), length.out = ncomp)

    d <- grid[2] - grid[1]

    g_init <- ashr::normalmix(pi = rep(1 / ncomp, ncomp),
                              mean = grid,
                              sd = d / 2)
  }

  return(g_init)
}

#' @importFrom stats approx
#'
default_smn_scale <- function(x,
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

#' @importFrom stats approx
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @importFrom dplyr filter arrange pull
#'
default_symmuni_scale <- function(x,
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
    scale <- default_smn_scale(
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
