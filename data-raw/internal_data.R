devtools::load_all(".")
library(tidyverse)


m <- exp(exp(seq(log(log(1.025)), log(log(4)), length.out = 100)))
KL <- numeric(length(m))

cat("Computing scale mixture of normals grid")
for (i in 1:length(m)) {
  cat(".")
  if (i %% 50 == 0) {
    cat("\n")
  }
  KL[i] <- ub.smnKLdiv(m[i])$KL.div
}
cat("\n")

smngrid <- tibble(KL = KL, m = m)


max.K <- 300
symmunigrid <- tibble()
for (log.KL in seq(-2, -2.75, by = -0.25)) {
  grid <- build.symm.grid(10^(log.KL), max.K = max.K, init.iter = 5, init.srch = c(1, 10))
  symmunigrid <- symmunigrid %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-3, -4.75, by = -0.25)) {
  grid <- build.symm.grid(10^(log.KL), max.K = max.K, init.iter = 10, init.srch = c(0.2, 2))
  symmunigrid <- symmunigrid %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-5, -5.75, by = -0.25)) {
  grid <- build.symm.grid(10^(log.KL), max.K = max.K, init.iter = 15, init.srch = c(0.1, 1))
  symmunigrid <- symmunigrid %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-6, -6.75, by = -0.25)) {
  grid <- build.symm.grid(10^(log.KL), max.K = max.K, init.iter = 25, init.srch = c(0.05, 1))
  symmunigrid <- symmunigrid %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}
for (log.KL in seq(-7, -8, by = -0.25)) {
  grid <- build.symm.grid(10^(log.KL), max.K = max.K, init.iter = 50, init.srch = c(0.01, 0.5))
  symmunigrid <- symmunigrid %>%
    bind_rows(tibble(KL = 10^(log.KL), idx = 1:length(grid), loc = grid))
}


usethis::use_data(smngrid, symmunigrid, internal = TRUE, overwrite = TRUE)
