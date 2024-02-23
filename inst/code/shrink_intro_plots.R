# TO DO: Explain what this script is for, and how to use it.
library(rmarkdown)
library(ggplot2)
library(cowplot)

# Run the shrink_intro vignette.
render("shrink_intro.Rmd")

# Compile the results needed to generate the plots.
pdat <- data.frame(u     = u,
                   mle   = x,
                   est.n = coef(fit_normal)
)

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

print(plot_grid(p1,p2,nrow = 1,ncol = 2))
