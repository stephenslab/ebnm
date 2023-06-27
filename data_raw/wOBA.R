library(tidyverse)

yr <- 2022

# Get wOBA weights for season of interest:
fg_guts <- baseballr::fg_guts()
w <- fg_guts %>%
  filter(season == yr) %>%
  select(wBB:wHR) %>%
  as.matrix()

# Select relevant variables:
fg_data_all <- baseballr::fg_batter_leaders(
  yr, yr, qual = 1, ind = 1, exc_p = TRUE
)
fg_data <- fg_data_all %>%
  select(playerid, Name, Team, PA, BB, HBP, `1B`, `2B`, `3B`, HR, wOBA)

# Calculate the wOBA variance for 1 PA that would be obtained using league-average
#.  outcome proportions:
lg_props <- fg_data %>%
  summarize(across(PA:HR, sum)) %>%
  mutate(across(BB:HR, ~ . / PA)) %>%
  select(-PA) %>%
  as.matrix()
lg_covmat <- - t(lg_props) %*% lg_props
diag(lg_covmat) <- lg_props * (1 - lg_props)
lg_wOBA_var <- w %*% lg_covmat %*% t(w)

# Calculate wOBA variance for 1 PA using empirical outcome proportions:
estmat <- fg_data %>%
  mutate(across(BB:HR, ~ . / PA)) %>%
  select(BB:HR) %>%
  as.matrix()
plugin_wOBA_var <- sapply(1:nrow(fg_data), function(i) {
  pihat <- estmat[i, , drop = FALSE]
  covmat <- -t(pihat) %*% pihat
  diag(covmat) <- pihat * (1 - pihat)
  return(w %*% covmat %*% t(w))
})

# For each hitter, use the greater of the two estimates:
wOBA_var <- pmax(plugin_wOBA_var, lg_wOBA_var)

wOBA <- fg_data %>%
  transmute(
    FanGraphsID = playerid,
    Name,
    Team,
    PA = PA,
    x = wOBA,
    s = sqrt(wOBA_var / PA)
  )

usethis::use_data(wOBA, overwrite = TRUE)
