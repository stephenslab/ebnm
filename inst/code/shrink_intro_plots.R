# This script generates the plots for Fig. 1 of the paper. This should
# be run in the vignettes directory of the ebnm package source, e.g.,
#
#   source("../inst/code/shrink_intro_plots.R")
#
library(rmarkdown)
library(ggplot2)
library(cowplot)

# Run the shrink_intro vignette.
render("shrink_intro.Rmd")

# Compile the results needed to generate the plots.
pdat <- data.frame(u      = u,
                   mle    = x,
                   est.n  = coef(fit_normal),
                   est.pn = coef(fit_pn))

# Plot true mean vs. MLE.
lims <- c(-0.55,5.05)
p1 <- ggplot(pdat,aes(x = u,y = mle)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  labs(x = "true value",
       y = "estimate",
       title = "MLE",
       subtitle = sprintf("MSE = %0.3f",mean(err_mle))) +
  xlim(lims) +
  ylim(lims) +
  theme_cowplot(font_size = 12) +
  theme(plot.title = element_text(face = "plain",size = 12),
        plot.subtitle = element_text(face = "plain",size = 12))

# Plot true mean vs. posterior estimate from fit_normal.
p2 <- ggplot(pdat,aes(x = u,y = est.n)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  labs(x = "true value ",
       y = "estimate",
       title = "EB with normal prior",
       subtitle = sprintf("MSE = %0.3f",mean(err_shrink_normal))) +
  xlim(lims) +
  ylim(lims) +
  theme_cowplot(font_size = 12) +
  theme(plot.title = element_text(face = "plain",size = 12),
        plot.subtitle = element_text(face = "plain",size = 12))

# Plot true mean vs. posterior estimate from fit_pn.
p3 <- ggplot(pdat,aes(x = u,y = est.pn)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  labs(x = "true value ",
       y = "estimate",
       title = "EB with point-normal prior",
       subtitle = sprintf("MSE = %0.3f",mean(err_shrink_pn))) +
  xlim(lims) +
  ylim(lims) +
  theme_cowplot(font_size = 12) +
  theme(plot.title = element_text(face = "plain",size = 12),
        plot.subtitle = element_text(face = "plain",size = 12))

print(plot_grid(p1,p2,p3,nrow = 1,ncol = 3))
ggsave("shrink_intro.eps",
       plot_grid(p1,p2,p3,nrow = 1,ncol = 3),
       height = 2.75,width = 7)
