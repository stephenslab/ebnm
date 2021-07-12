## Quality of ashr grid approximations for scale mixtures of normals.

# To sample from different distributions N(0, s2), use a single sample from
#   N(0, 1) and scale. This method ensures a smooth optimization objective.
smnKLdiv <- function(s2, omega, m, samp) {
  samp <- samp * sqrt(s2)

  tmp <- omega * exp(0.5 * samp^2 * (1 / s2 - 1))
  tmp <- tmp + (1 - omega) * exp(0.5 * samp^2 * (1 / s2 - 1 / m)) / sqrt(m)
  tmp <- -log(tmp) - 0.5 * log(s2)

  return(mean(tmp))
}

min_smnKLdiv <- function(s2, m, samp) {
  optres <- optimize(
    function(omega) smnKLdiv(s2, omega, m, samp),
    interval = c(0, 1),
    maximum = FALSE
  )
  return(optres$objective)
}

ub_smnKLdiv <- function(m, samp) {
  optres <- optimize(
    function(s2) min_smnKLdiv(s2, m, samp),
    interval = c(1, m),
    maximum = TRUE
  )
  return(optres$objective)
}

build_grid <- function(m, sampsize) {
  cat("Building grid for values from", min(m), "to", max(m), "\n")
  samp <- rnorm(sampsize)

  ub <- numeric(length(m))
  for (i in 1:length(m)) {
    cat("  m:", m[i], "\n")
    ub[i] <- ub_smnKLdiv(m[i], samp)
  }

  tib <- tibble::tibble(m = m, ub = ub)

  return(tib)
}

set.seed(666)

# The relationship between log(ub) and log(log(m)) is roughly linear.
tib <- build_grid(
  m = exp(exp(seq(log(log(1.05)), log(log(4)), length.out = 100))),
  sampsize = 1e8
)

saveRDS(tib, "./data/smngrid.rds")
