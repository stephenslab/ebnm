#' @importFrom stats pnorm
#'
pl_nlm_fn <- function(par, x, s, lf) {
  w <- 1 / (1 + exp(-par[1]))
  a <- exp(par[2])

  # This part of the log likelihood comes from the left tail.
  xleft <- x / s + s * a
  lpnormleft <- pnorm(xleft, log.p = TRUE, lower.tail = FALSE)
  lgleft <- a * x + lpnormleft

  # This part comes from the right tail.
  xright <- x / s - s * a
  lpnormright <- pnorm(xright, log.p = TRUE)
  lgright <- -a * x + lpnormright

  # This part corresponds to the entire Laplace component.
  lg <- log(a / 2) + s^2 * a^2 / 2 + logscale_add(lgleft, lgright)

  llik <- logscale_add(log(1 - w) + lf, log(w) + lg)

  # Derivatives with respect to w and alpha (= par[1]).
  f.over.lik <- exp(lf - llik)
  g.over.lik <- exp(lg - llik)
  dnllik.dw <- f.over.lik - g.over.lik
  dnllik.dalpha <- w * (1 - w) * dnllik.dw
  d2nllik.dalpha2 <- (1 - 2 * w) * dnllik.dalpha + dnllik.dalpha^2

  # Derivatives with respect to a and beta (= par[2]).
  dlpnormleft.da <- -s * exp(-log(2 * pi) / 2 - xleft^2 / 2 - lpnormleft)
  dlgleft.da <- x + dlpnormleft.da
  d2lgleft.da2 <- dlpnormleft.da * (-s * xleft - dlpnormleft.da)

  dlpnormright.da <- -s * exp(-log(2 * pi) / 2 - xright^2 / 2 - lpnormright)
  dlgright.da <- -x + dlpnormright.da
  d2lgright.da2 <- dlpnormright.da * (s * xright - dlpnormright.da)

  wt <- 1 / (1 + exp(lgleft - lgright))
  dlg.da <- 1 / a + s^2 * a + (1 - wt) * dlgleft.da + wt * dlgright.da
  d2lg.da2 <- -1 / a^2 + s^2 + (1 - wt) * d2lgleft.da2 + wt * d2lgright.da2
  d2lg.da2 <- d2lg.da2 + wt * (1 - wt) * (dlgleft.da - dlgright.da)^2

  dnllik.da <- -w * g.over.lik * dlg.da
  d2nllik.da2 <- -w * g.over.lik * (d2lg.da2 + dlg.da * (dlg.da + dnllik.da))

  dnllik.dbeta <- a * dnllik.da
  d2nllik.dbeta2 <- a^2 * d2nllik.da2 + dnllik.dbeta

  # Mixed second derivative.
  d2nllik.dalphadbeta <- (dnllik.dalpha * dnllik.dbeta
                          - (w * (1 - w) * a) * g.over.lik * dlg.da)

  retval <- -sum(llik)
  attr(retval, "gradient") <- c(sum(dnllik.dalpha), sum(dnllik.dbeta))
  attr(retval, "hessian") <- matrix(c(sum(d2nllik.dalpha2),
                                      rep(sum(d2nllik.dalphadbeta), 2),
                                      sum(d2nllik.dbeta2)),
                                    nrow = 2, ncol = 2)
  return(retval)
}

logscale_add <- function(log.x, log.y) {
  C <- pmax(log.x, log.y)
  return(log(exp(log.x - C) + exp(log.y - C)) + C)
}
