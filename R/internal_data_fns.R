## Scale mixture of normal grids ----------------------------------------------

log_add <- function(log.x, log.y) {
  c <- pmax(log.x, log.y)
  return(log(exp(log.x - c) + exp(log.y - c)) + c)
}

# KL divergence from f = N(0, s2) to g = omega * N(0, 1) + (1 - omega) * N(0, m).
#
#' @importFrom stats dnorm integrate
#'
smnKLdiv <- function(s2, omega, m) {
  KL.integrand <- function(x) {
    f.dens <- dnorm(x, mean = 0, sd = sqrt(s2))
    logf.dens <- dnorm(x, mean = 0, sd = sqrt(s2), log = TRUE)
    logg1.dens <- log(omega) + dnorm(x, mean = 0, sd = 1, log = TRUE)
    logg2.dens <- log(1 - omega) + dnorm(x, mean = 0, sd = sqrt(m), log = TRUE)
    logg.dens <- log_add(logg1.dens, logg2.dens)
    return(f.dens * (logf.dens - logg.dens))
  }
  int.res <- integrate(
    KL.integrand,
    lower = -Inf,
    upper = Inf,
    rel.tol = sqrt(.Machine$double.eps)
  )
  return(int.res$value)
}

# Minimum KL divergence over 0 \le omega \le 1.
#
#' @importFrom stats optimize
#'
min_smnKLdiv <- function(s2, m) {
  optres <- optimize(
    function(omega) smnKLdiv(s2, omega, m),
    interval = c(0, 1),
    maximum = FALSE
  )
  return(optres)
}

# Maximum KL divergence over 1 \le s2 \le m.
#
#' @importFrom stats optimize
#'
ub.smnKLdiv <- function(m) {
  optres <- optimize(
    function(s2) min_smnKLdiv(s2, m)$objective,
    interval = c(1, m),
    maximum = TRUE
  )
  retlist <- list(
    KL.div = optres$objective,
    opt.s2 = optres$maximum,
    opt.omega = min_smnKLdiv(optres$maximum, m)$minimum
  )
  return(retlist)
}


## Symmetric unimodal grids ---------------------------------------------------

# Returns log(x - y), so assumes x > y.
log_minus <- function(log.x, log.y) {
  return(log(1 - exp(log.y - log.x)) + log.x)
}

# Log density of UN(a, 1) (Unif[-a, a] convolved with N(0, 1)).
#
#' @importFrom stats dnorm pnorm
#'
UN.logdens <- function(x, a) {
  x <- abs(x) # computations are stabler for x > 0
  if (a == 0) {
    retval <- dnorm(x, log = TRUE)
  } else {
    retval <- log_minus(
      pnorm(a - x, log.p = TRUE),
      pnorm(-a - x, log.p = TRUE)
    ) - log(2 * a)
  }
  return(retval)
}

# KL divergence from f = UN(a, 1) to
#   g = omega * UN(a_left, 1) + (1 - omega) * UN(a_right, 1).
#
#' @importFrom stats integrate
#'
symmKLdiv <- function(a, omega, a.left, a.right) {
  KL.integrand <- function(x) {
    f.dens <- exp(UN.logdens(x, a = a))
    logf.dens <- UN.logdens(x, a = a)
    logg1.dens <- log(omega) + UN.logdens(x, a = a.left)
    logg2.dens <- log(1 - omega) + UN.logdens(x, a = a.right)
    logg.dens <- log_add(logg1.dens, logg2.dens)
    return(f.dens * (logf.dens - logg.dens))
  }
  int.res <- integrate(
    KL.integrand,
    lower = -a - 6,
    upper = a + 6,
    rel.tol = sqrt(.Machine$double.eps)
  )
  return(int.res$value)
}

# Minimum KL divergence over 0 \le omega \le 1.
#
#' @importFrom stats optimize
#'
min_symmKLdiv <- function(a, a.left, a.right) {
  optres <- optimize(
    function(omega) symmKLdiv(a, omega, a.left, a.right),
    interval = c(0, 1),
    maximum = FALSE
  )
  return(optres)
}

# Maximum KL divergence over a_left \le a \le a_right.
#
#' @importFrom stats optimize
#'
ub.symmKLdiv <- function(a.left, a.right) {
  optres <- optimize(
    function(a) min_symmKLdiv(a, a.left, a.right)$objective,
    interval = c(a.left, a.right),
    maximum = TRUE
  )
  retlist <- list(
    KL.div = optres$objective,
    opt.a = optres$maximum,
    opt.omega = min_symmKLdiv(optres$maximum, a.left, a.right)$minimum
  )
  return(retlist)
}

# Given a_left and a desired (upper bound for) KL-divergence, finds a_right.
#
#' @importFrom stats uniroot
#'
symm.find.next.gridpt <- function(a.left, targetKL, srch.interval) {
  uniroot.fn <- function(a.right) {
    return(ub.symmKLdiv(a.left, a.right)$KL.div - targetKL)
  }
  optres <- uniroot(
    uniroot.fn,
    c(a.left + srch.interval[1], a.left + srch.interval[2])
  )
  return(optres$root)
}

# Builds a grid starting from a_left = 0. Uses init.srch.interval for the first
#   init.iter iterations, then assumes that a_right / a_left decreases to its theoretical limit.
build.symm.grid <- function(targetKL,
                       max.K = 30,
                       max.x = Inf,
                       init.srch.interval = c(0.1, 5),
                       init.iter = 5) {
  print.progress <- function(K) {
    if (K %% 50 == 0) {
      cat("X\n")
    } else if (K %% 10 == 0) {
      cat("X")
    } else {
      cat(".")
    }
  }

  cat("Building grid ( KL =", targetKL, ")..")
  grid <- 0
  K <- 1

  is.space.increasing <- FALSE
  while (K < init.iter && max(grid) < max.x) {
    next.gridpt <- symm.find.next.gridpt(
      grid[K],
      targetKL,
      srch.interval = init.srch.interval
    )
    grid <- c(grid, next.gridpt)
    K <- K + 1
    print.progress(K)
  }

  min.incr <- 1 + exp(1) * targetKL
  while (length(grid) < max.K && max(grid) < max.x) {
    last.space <- grid[length(grid)] - grid[length(grid) - 1]
    last.incr <- grid[length(grid)] / grid[length(grid) - 1]
    next.gridpt <- symm.find.next.gridpt(
      grid[K],
      targetKL,
      srch.interval = c(last.space * min.incr, last.space * last.incr)
    )
    grid <- c(grid, next.gridpt)
    K <- K + 1
    print.progress(K)
  }
  cat("\n")
  return(grid)
}


## Grids for NPMLE ------------------------------------------------------------

# These functions aren't used because good, simple upper bounds exist.

# KL divergence from f = N(0, 1) to g = 0.5 * N(-d/2, s2 + 1) + 0.5 * N(d/2, s2 + 1):
#
#' @importFrom stats dnorm integrate
#'
npmleKLdiv <- function(d, s2) {
  KL.integrand <- function(x) {
    f.dens <- dnorm(x, mean = 0, sd = 1)
    logf.dens <- dnorm(x, mean = 0, sd = 1, log = TRUE)
    logg1.dens <- log(0.5) + dnorm(x, mean = -d / 2, sd = sqrt(s2 + 1), log = TRUE)
    logg2.dens <- log(0.5) + dnorm(x, mean = d / 2, sd = sqrt(s2 + 1), log = TRUE)
    logg.dens <- log_add(logg1.dens, logg2.dens)
    return(f.dens * (logf.dens - logg.dens))
  }
  int.res <- integrate(
    KL.integrand,
    lower = -Inf,
    upper = Inf,
    rel.tol = sqrt(.Machine$double.eps)
  )
  return(int.res$value)
}

# Minimum KL divergence over s2 \ge 0.
#
#' @importFrom stats optimize
#'
min_npmleKLdiv <- function(d) {
  optres <- optimize(
    function(s2) npmleKLdiv(d, s2),
    interval = c(0, d^2),
    maximum = FALSE
  )
  return(optres)
}
