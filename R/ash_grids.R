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
                              mode,
                              min_K = 3,
                              max_K = 300,
                              KLdiv_target = 1 / length(x)) {
  max_x2 <- max((x - mode)^2)
  min_s2 <- min(s)^2

  max_mult <- (max(max_x2, min_s2) + 1)^(1 / (min_K - 1))
  min_mult <- max_mult^((min_K - 1) / (max_K - 1))

  grid_mult <- approx(
    y = smngrid$m,
    x = log(smngrid$KL),
    xout = log(KLdiv_target),
    rule = 2
  )[["y"]]

  grid_mult <- min(max(grid_mult, min_mult), max_mult)

  K <- ceiling(log(max(max_x2 / min_s2, 1), base = grid_mult) + 1)
  scale <- sqrt(min_s2 * (grid_mult^(0:(K - 1)) - 1))

  return(scale)
}

#' @importFrom stats approx
#'
default_symmuni_scale <- function(x,
                                  s,
                                  mode,
                                  min_K = 3,
                                  max_K = 300,
                                  KLdiv_target = 1 / length(x)) {
  return(default_smn_scale(x, s, mode))
}

get_ashr_grid <- function(x, s, mode, grid_mult) {
  if (grid_mult == "default") {
    grid_mult <- sqrt(2)
  }

  # Adapted from ashr:::autoselect.mixsd.
  sigmamin <- min(s[s > 0]) / 10
  sigmamax <- max(8 * sigmamin, 2 * sqrt(max((x - mode)^2 - s^2, 0)))
  npoint <- ceiling(log2(sigmamax / sigmamin) / log2(grid_mult))
  scale <- grid_mult^((-npoint):0) * sigmamax
  scale <- c(0, scale)

  return(scale)
}
