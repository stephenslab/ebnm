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
  } else if (force_pointmass || xrange / max_K < 2.3 * min(s)) {
    # Use point-mass mixture.
    d <- min(s) * (64 * KLdiv_target)^0.25
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
default_scale <- function(x,
                          s,
                          mode,
                          min_K = 3,
                          max_K = 300,
                          KLdiv_target = 1 / length(x)) {
  max_x2 <- max((x - mode)^2)
  min_s2 <- min(s)^2

  max_mult <- (max(max_x2, min_s2) + 1)^(1 / (min_K - 1))
  min_mult <- max_mult^((min_K - 1) / (max_K - 1))

  # To get arguments for approx, run:
  #   readRDS("./data/smngrid.rds")
  #   paste0("y = c(", paste(signif(smngrid$m, digits = 2), collapse = ", "), ")")
  #   paste0("x = c(", paste(signif(log(smngrid$ub), digits = 2), collapse = ", "), ")")
  grid_mult <- approx(
    y = c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0),
    x = c(-14.0, -10.0, -9.2, -8.4, -7.8, -7.3, -7.0, -6.6, -6.4),
    xout = log(KLdiv_target),
    rule = 2
  )[["y"]]

  grid_mult <- min(max(grid_mult, min_mult), max_mult)

  K <- ceiling(log(max(max_x2 / min_s2, 1), base = grid_mult) + 1)
  K <- min(max(K, min_K), max_K)
  scale <- sqrt(min_s2 * (grid_mult^(0:(K - 1)) - 1))

  return(scale)
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
