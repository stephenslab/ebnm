library(ggplot2)
library(cowplot)
library(reshape2)
breaks <- seq(-0.6,4.5,length.out = 20)
dat <- data.frame(u,
                  mle = (x - u)^2,
                  shrink = (u - coef(fit_normal))^2,
                  shrink_pn = (u - coef(fit_pn))^2)
dat$u <- cut(dat$u,breaks)
pdat <- data.frame(u = (breaks[-1] + breaks[-16])/2,
                   mle = tapply(dat$mle,dat$u,mean),
                   shrink = tapply(dat$shrink,dat$u,mean),
                   shrink_pn = tapply(dat$shrink_pn,dat$u,mean))
pdat <- melt(pdat,id = "u",variable.name = "method")
p <- ggplot(pdat,aes(x = u,y = value,color = method)) +
  geom_line(size = 0.75) +
  labs(y = "mse") +
  theme_cowplot()
print(p)
