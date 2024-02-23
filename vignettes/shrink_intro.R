## ----knitr-opts, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment = "#",collapse = TRUE,results = "hold",
                      fig.align = "center",dpi = 90)


## ----sim-data-----------------------------------------------------------------
set.seed(1)
n <- 100000
u <- 2 + (runif(n) < 0.2) * rnorm(n)


## ----sim-data-2---------------------------------------------------------------
s <- rep(1/3,n)
x <- u + s*rnorm(n)


## ----plot-mle, fig.height=3.5, fig.width=3.5----------------------------------
par(mar = c(4,4,2,2))
lims <- range(c(u,x))
plot(u,x,pch = 4,cex = 0.75,xlim = lims,ylim = lims)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")


## ----ebnm-normal--------------------------------------------------------------
library(ebnm)
fit_normal <- ebnm(x,s,prior_family = "normal",mode = "estimate")


## ----plot-ebnm-normal, fig.height=3.5, fig.width=3.5--------------------------
y <- coef(fit_normal)
lims <- range(c(u,y))
par(mar = c(4,4,2,2))
plot(u,y,pch = 4,cex = 0.75,xlim = lims,ylim = lims)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")


## ----mse-1--------------------------------------------------------------------
mse_mle    <- (x - u)^2
mse_shrink <- (y - u)^2
print(round(digits = 4,
            x = c(mle = mean(mse_mle),
                  normal = mean(mse_shrink))))


## ----plot-mse-1, fig.height=3.5, fig.width=3.5--------------------------------
par(mar = c(4,4,2,2))
plot(mse_mle,mse_shrink,pch = 4,cex = 0.75,xlim = c(0,1.2),ylim = c(0,1.2))
abline(a = 0,b = 1,col = "magenta",lty = "dotted")


## ----ebnm-pn------------------------------------------------------------------
fit_pn <- ebnm(x,s,prior_family = "point_normal",mode = "estimate")


## ----plot-ebnm-pn, fig.height=3.5, fig.width=3.5------------------------------
par(mar = c(4,4,2,2))
y <- coef(fit_pn)
lim <- range(c(u,y))
plot(u,y,pch = 4,cex = 0.75,xlim = lim,ylim = lim)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")


## ----plot-mse-2, fig.height=3.5, fig.width=3.5--------------------------------
mse_shrink_pn <- (y - u)^2
print(round(digits = 4,
            x = c(mle = mean(mse_mle),
                  normal = mean(mse_shrink),
				  point_normal = mean(mse_shrink_pn))))

