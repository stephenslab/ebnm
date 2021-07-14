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
  #   smngrid <- readRDS("./data/smngrid.rds")
  #   paste0("y = c(", paste(signif(smngrid$m, digits = 6), collapse = ", "), "),")
  #   paste0("x = c(", paste(signif(log(smngrid$ub), digits = 6), collapse = ", "), "),")
  grid_mult <- approx(
    y = c(1.05, 1.05176, 1.05359, 1.05548, 1.05744, 1.05948, 1.06158, 1.06377, 1.06603, 1.06838, 1.07081, 1.07333, 1.07595, 1.07866, 1.08147, 1.08439, 1.08741, 1.09055, 1.0938, 1.09718, 1.10069, 1.10432, 1.1081, 1.11202, 1.11608, 1.12031, 1.12469, 1.12924, 1.13397, 1.13889, 1.14399, 1.14929, 1.15481, 1.16054, 1.16649, 1.17268, 1.17913, 1.18582, 1.19279, 1.20005, 1.2076, 1.21545, 1.22364, 1.23216, 1.24103, 1.25028, 1.25992, 1.26997, 1.28045, 1.29138, 1.30279, 1.31469, 1.32712, 1.34009, 1.35365, 1.36782, 1.38263, 1.39812, 1.41432, 1.43128, 1.44904, 1.46763, 1.48712, 1.50755, 1.52898, 1.55147, 1.57508, 1.59988, 1.62594, 1.65334, 1.68217, 1.71253, 1.7445, 1.7782, 1.81374, 1.85126, 1.89088, 1.93276, 1.97705, 2.02393, 2.07359, 2.12625, 2.18212, 2.24146, 2.30454, 2.37166, 2.44314, 2.51934, 2.60067, 2.68756, 2.78049, 2.88, 2.98667, 3.10118, 3.22425, 3.35668, 3.4994, 3.65341, 3.81985, 4),
    x = c(-16.1427, -16.0288, -15.9139, -15.7982, -15.6816, -15.5641, -15.4458, -15.3267, -15.2068, -15.0862, -14.9648, -14.8428, -14.7201, -14.5969, -14.473, -14.3485, -14.2236, -14.0981, -13.9721, -13.8457, -13.7189, -13.5917, -13.4641, -13.3362, -13.208, -13.0794, -12.9506, -12.8216, -12.6923, -12.5629, -12.4332, -12.3034, -12.1735, -12.0434, -11.9132, -11.783, -11.6527, -11.5223, -11.392, -11.2616, -11.1312, -11.0009, -10.8706, -10.7403, -10.6102, -10.4801, -10.3501, -10.2203, -10.0906, -9.9611, -9.83176, -9.70261, -9.57367, -9.44496, -9.31649, -9.18829, -9.06035, -8.93271, -8.80538, -8.67837, -8.55171, -8.4254, -8.29946, -8.17391, -8.04876, -7.92403, -7.79974, -7.67589, -7.55251, -7.42961, -7.3072, -7.1853, -7.06391, -6.94306, -6.82276, -6.70301, -6.58384, -6.46525, -6.34725, -6.22986, -6.11308, -5.99693, -5.88142, -5.76655, -5.65233, -5.53877, -5.42588, -5.31367, -5.20214, -5.0913, -4.98116, -4.87172, -4.76298, -4.65496, -4.54766, -4.44107, -4.33521, -4.23009, -4.12569, -4.02203),
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
