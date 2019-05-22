likelihood.include <- 'Rcpp::NumericVector llik(SEXP xs, SEXP env) {
  Rcpp::NumericVector par(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);

  int fixpi0 = e["fixpi0"];
  int fixa = e["fixa"];
  int fixmu = e["fixmu"];
  int i = 0;

  double n0 = e["n0"];
  double n1 = e["n1"];
  double n2 = e["n2"];
  arma::vec s2 = Rcpp::as<arma::vec>(e["s2"]);

  double alpha;
  if (fixpi0) {
    alpha = e["alpha"];
  } else {
    alpha = par[i]; i++;
  }

  double beta;
  if (fixa) {
    beta = e["beta"];
  } else {
    beta = par[i]; i++;
  }

  double mu;
  double sum1;
  arma::vec z;
  if (fixmu) {
    mu = e["mu"];
    sum1 = e["sum1"];
    z = Rcpp::as<arma::vec>(e["z"]);
  } else {
    mu = par[i];
    sum1 = sum(square(Rcpp::as<arma::vec>(e["x1"]) - mu));
    z = square(Rcpp::as<arma::vec>(e["x2"]) - mu) / s2;
  }

  double nllik0 = n0 * log(1 + exp(-alpha));
  double nllik1 = n1 * (log(1 + exp(alpha))) + n1 * beta / 2 + sum1 * exp(-beta) / 2;
  arma::vec y = (z / (1 + s2 * exp(-beta)) - log(1 + exp(beta) / s2)) / 2;
  arma::vec C = max(y, alpha * arma::ones<arma::vec>(n2));
  double nllik2 = n2 * (log(1 + exp(alpha))) - sum(log(exp(y - C) + exp(alpha - C)) + C);
  double out = nllik0 + nllik1 + nllik2;

  Rcpp::NumericVector ret = Rcpp::as<Rcpp::NumericVector>(wrap(out));
  return ret;
}'

gradient.include <- 'Rcpp::NumericVector grad(SEXP xs, SEXP env) {
  Rcpp::NumericVector par(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);

  int fixpi0 = e["fixpi0"];
  int fixa = e["fixa"];
  int fixmu = e["fixmu"];
  int i = 0;

  double n0 = e["n0"];
  double n1 = e["n1"];
  double n2 = e["n2"];
  arma::vec s2 = Rcpp::as<arma::vec>(e["s2"]);

  double alpha;
  if (fixpi0) {
    alpha = e["alpha"];
  } else {
    alpha = par[i]; i++;
  }

  double beta;
  if (fixa) {
    beta = e["beta"];
  } else {
    beta = par[i]; i++;
  }

  double mu;
  double sum1;
  arma::vec x;
  arma::vec s;
  arma::vec z;
  arma::vec x1;
  if (fixmu) {
    mu = e["mu"];
    sum1 = e["sum1"];
    z = Rcpp::as<arma::vec>(e["z"]);
  } else {
    x = Rcpp::as<arma::vec>(e["x"]);
    s = Rcpp::as<arma::vec>(e["s"]);
    mu = par[i];
    x1 = Rcpp::as<arma::vec>(e["x1"]);
    sum1 = sum(square(x1 - mu));
    z = square(x - mu) / s2;
  }

  arma::vec tmp1 = 1 / (1 + s2 * exp(-beta));
  arma::vec tmp2 = 1 / (1 + exp(beta) / s2);
  arma::vec y = (z % tmp1 + log(tmp2)) / 2;

  arma::vec grad(3 - fixpi0 - fixa - fixmu);

  int j = 0;
  double tmp;
  if (!fixpi0) {
    tmp = -n0 / (1 + exp(alpha)) + (n1 + n2) / (1 + exp(-alpha));
    grad[j] = tmp - sum(1 / (1 + exp(y - alpha)));
    j++;
  }

  if (!fixa) {
    tmp = n1 / 2 - exp(-beta) * sum1 / 2;
    grad[j] = tmp - sum((z % tmp1 % tmp2 - tmp1) / (1 + exp(alpha - y))) / 2;
    j++;
  }

  if (!fixmu) {
    tmp = sum(mu - x1) * exp(-beta);
    grad[j] = tmp + sum(tmp1 % (x - mu) / (1 + exp(alpha - y)) / s);
  }

  Rcpp::NumericVector ret = Rcpp::as<Rcpp::NumericVector>(wrap(grad));
  return ret;
}'

likelihood.body <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&llik)));
'

gradient.body <- '
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&grad)));
'

likelihood.CPP <- inline::cxxfunction(signature(), body=likelihood.body,
                                      inc=likelihood.include, plugin="RcppArmadillo")

gradient.CPP <- inline::cxxfunction(signature(), body=gradient.body,
                                    inc=gradient.include, plugin="RcppArmadillo")

cpp_mle_point_normal <- function(x, s, g, control, fix_pi0, fix_mu) {
  startpar <- rep(0, 3 - fix_pi0 - fix_mu)
  i <- 1
  if (!is.null(g)) {
    if (!is.null(g$pi0) && !fix_pi0) {
      startpar[i] <- -log(1 / g$pi0 - 1)
      i <- i + 1
    }
    if (!is.null(g$a)) {
      startpar[i] <- -log(g$a)
      i <- i + 1
    }
    if (!is.null(g$mu) && !fix_mu) {
      startpar[i] <- g$mu
    }
  }

  env <- new.env()
  env[["fixpi0"]] <- 1L * fix_pi0
  env[["fixa"]] <- 0L
  env[["fixmu"]] <- 1L * fix_mu
  if (fix_pi0) {
    env[["alpha"]] <- -log(1 / g$pi0 - 1)
  }
  if (fix_mu) {
    env[["mu"]] <- g$mu
  }
  if (any(s == 0)) {
    which.s0 <- which(s == 0)
    s <- s[-which.s0]
    if (fix_mu) {
      which.x.nz <- which(x[which.s0] != 0)
      env[["n0"]] <- length(which.s0) - length(which.x.nz)
      env[["n1"]] <- length(which.x.nz)
      env[["sum1"]] <- sum(x[which.x.nz]^2)
    } else {
      env[["n0"]] <- 0
      env[["n1"]] <- length(which.s0)
      env[["sum1"]] <- sum(x[which.s0]^2)
      env[["x1"]] <- x[which.s0]
      env[["s"]] <- s
    }
    x <- x[-which.s0]
    env[["x"]] <- x
 } else {
    env[["n0"]] <- 0
    env[["n1"]] <- 0
    env[["sum1"]] <- 0
  }
  env[["n2"]] <- length(x)
  if (length(s) == 1) {
    env[["s2"]] <- rep(s^2, length(x))
  } else {
    env[["s2"]] <- s^2
  }
  env[["z"]] <- (x / s)^2

  opt.res <- lbfgs::lbfgs(likelihood.CPP(), gradient.CPP(),
                          startpar, environment = env, invisible = 1)
  out <- list()
  i <- 1
  if (fix_pi0) {
    out$pi0 <- g$pi0
  } else {
    out$pi0 <- 1 / (exp(-opt.res$par[i]) + 1)
    i <- i + 1
  }
  out$a <- exp(-opt.res$par[i])
  i <- i + 1
  if (fix_mu) {
    out$mu <- g$mu
  } else {
    out$mu <- opt.res$par[i]
  }

  out$val <- opt.res$value + 0.5 * ((env[["n1"]] + env[["n2"]]) * log(2 * pi)
                                    + sum(log(env[["s2"]])) + sum(env[["z"]]))

  return(out)
}


# Functions to compute MLE under point-normal prior. Essentially a wrapper
#   to optim.

mle_point_normal_logscale_grad <- function(x, s, g, control, fix_pi0, fix_mu) {
  # Optimization functions. "par" is the set of variables that is being
  #   optimized: -logit(pi0) (if fix_pi0 is FALSE), log(a), and mu (if fix_mu
  #   is FALSE).
  fn <- function(par) {
    return(-loglik_point_normal(x,
                                s,
                                w = w_from_par(par, g, fix_pi0),
                                a = a_from_par(par, fix_pi0),
                                mu = mu_from_par(par, g, fix_pi0, fix_mu)))
  }
  gr <- function(par) {
    ret <- grad_negloglik_logscale_point_normal(x,
                                                s,
                                                w = w_from_par(par, g, fix_pi0),
                                                a = a_from_par(par, fix_pi0),
                                                mu = mu_from_par(par, g, fix_pi0, fix_mu))
    # Remove variables that aren't being optimized.
    return(ret[c(!fix_pi0, TRUE, !fix_mu)])
  }

  # Initial values.
  startpar <- numeric(0)
  if (!fix_pi0) {
    if (!is.null(g$pi0) && g$pi0 > 0 && g$pi0 < 1) {
      startpar <- c(startpar, log(1 / g$pi0 - 1))
    } else {
      startpar <- c(startpar, 0) # default for -logit(pi0)
    }
  }
  if (!is.null(g$a)) {
    startpar <- c(startpar, log(g$a))
  } else {
    startpar <- c(startpar, -log(mean(x^2))) # default for log(a)
  }
  if (!fix_mu) {
    if (!is.null(g$mu)) {
      startpar <- c(startpar, g$mu)
    } else {
      startpar <- c(startpar, mean(x)) # default for mu
    }
  }

  uu <- optimize_it(startpar, fn, gr, control,
                    mle_point_normal_hilo(x, s, fix_pi0, fix_mu),
                    x, s)

  # Pull pi0, a, and mu out of the optim results.
  retlist <- list()
  if (fix_pi0) {
    retlist$pi0 <- g$pi0
    retlist$a <- exp(uu$par[1])
    if (fix_mu) {
      retlist$mu <- g$mu
    } else {
      retlist$mu <- uu$par[2]
    }
  } else {
    retlist$pi0 <- 1 / (1 + exp(uu$par[1]))
    retlist$a <- exp(uu$par[2])
    if (fix_mu) {
      retlist$mu <- g$mu
    } else {
      retlist$mu <- uu$par[3]
    }
  }
  retlist$val <- uu$value
  return(retlist)
}

w_from_par <- function(par, g, fix_pi0) {
  if (fix_pi0)
    return(1 - g$pi0)
  else
    return(1 - (1/(1 + exp(par[1]))))
}

a_from_par <- function(par, fix_pi0) {
  if (fix_pi0)
    return(exp(par[1]))
  else
    return(exp(par[2]))
}

mu_from_par <- function(par, g, fix_pi0, fix_mu) {
  if (fix_pi0 && !fix_mu)
    return(par[2])
  else if (!fix_mu)
    return(par[3])
  else
    return(g$mu)
}


#' @importFrom stats optim
#'
optimize_it <- function(startpar, fn, gr, control, hilo, x, s) {
  uu <- try(lbfgsb3c::lbfgsb3c(startpar, fn, gr, control = control),
            silent=TRUE)

  # If optimization fails, try again with some limits; this should not
  # really happen but in preliminary testing sometimes we see optim
  # complain of infinite values, possibly because of extreme values of
  # the parameters?

  if (class(uu) == "try-error") {
    uu <- try(lbfgsb3c::lbfgsb3c(startpar, fn, gr,
                                 lower = hilo$lo, upper = hilo$hi, control = control))
  }

  # If optimization fails twice, save debug information and give up.

  if (class(uu) == "try-error") {
    saveRDS(list(startpar = startpar, x = x, s = s, control = control),
            "temp_debug.RDS")
    stop(paste("optim failed to converge; debug information saved to",
               "temp_debug.RDS"))
  }

  return(uu)
}


# Upper and lower bounds for optim in case the first attempt at optimization
#   fails.
mle_point_normal_hilo <- function(x, s, fix_pi0, fix_mu) {
  maxvar <- max(x^2)

  minvar <- (min(s) / 10)^2
  if (minvar < 1e-8) {
    minvar <- 1e-8
  }

  # Bounds for log(a):
  lo <- -log(maxvar)
  hi <- -log(minvar)

  if (!fix_pi0) {
    n <- length(x)
    lo <- c(log(1 / n), lo)
    hi <- c(log(n), hi)
  }

  if (!fix_mu) {
    lo <- c(lo, min(x) - 3 * max(s) - 3 * max(abs(x)))
    hi <- c(hi, max(x) + 3 * max(s) + 3 * max(abs(x)))
  }

  return(list(lo = lo, hi = hi))
}
