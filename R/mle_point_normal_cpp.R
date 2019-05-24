cpp_mle_point_normal <- function(x, s, g, control, fix_pi0, fix_a, fix_mu) {
  startpar <- pn_startpar(x, s, g, fix_pi0, fix_a, fix_mu)

  env <- new.env()
  env[["fixpi0"]] <- 1L * fix_pi0
  env[["fixa"]] <- 1L * fix_a
  env[["fixmu"]] <- 1L * fix_mu
  if (fix_pi0) {
    env[["alpha"]] <- -log(1 / g$pi0 - 1)
  }
  if (fix_a) {
    env[["beta"]] <- -log(g$a)
  }
  if (fix_mu) {
    env[["mu"]] <- g$mu
  }
  if (any(s == 0)) {
    which.s0 <- which(s == 0)
    which.x.nz <- which(x[which.s0] != 0)
    env[["n0"]] <- length(which.s0) - length(which.x.nz)
    env[["n1"]] <- length(which.x.nz)
    env[["sum1"]] <- sum(x[which.s0[which.x.nz]]^2)
    x <- x[-which.s0]
    s <- s[-which.s0]
  } else {
    env[["n0"]] <- 0
    env[["n1"]] <- 0
    env[["sum1"]] <- 0
  }
  env[["n2"]] <- length(x)
  env[["x"]] <- x
  if (length(s) == 1) {
    env[["s2"]] <- rep(s^2, length(x))
  } else {
    env[["s2"]] <- s^2
  }
  if (fix_mu) {
    env[["z"]] <- ((x - g$mu) / s)^2
  } else {
    env[["z"]] <- NULL
  }

  optres <- lbfgs::lbfgs(likelihood.CPP(),
                         gradient.CPP(),
                         startpar,
                         environment = env,
                         invisible = 1)

  retlist <- pn_g_from_optpar(optres$par, g, fix_pi0, fix_a, fix_mu)
  retlist$val <- pn_llik_from_optval(optres$value, env[["n1"]], env[["n2"]],
                                     env[["s2"]], env[["z"]])

  return(retlist)
}

likelihood.include <- 'Rcpp::NumericVector llik(SEXP xs, SEXP env) {
  Rcpp::NumericVector par(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);

  int fixpi0 = e["fixpi0"];
  int fixa = e["fixa"];
  int fixmu = e["fixmu"];

  double n0 = e["n0"];
  double n1 = e["n1"];
  double n2 = e["n2"];
  double sum1 = e["sum1"];
  arma::vec x = Rcpp::as<arma::vec>(e["x"]);
  arma::vec s2 = Rcpp::as<arma::vec>(e["s2"]);

  int i = 0;
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
  if (fixmu) {
    mu = e["mu"];
  } else {
    mu = par[i];
  }

  arma::vec z;
  if (fixmu) {
    z = Rcpp::as<arma::vec>(e["z"]);
  } else {
    z = square(x - mu) / s2;
  }
  arma::vec y = (z / (1 + s2 * exp(-beta)) - log(1 + exp(beta) / s2)) / 2;
  arma::vec C = max(y, alpha * arma::ones<arma::vec>(n2));

  double nllik = n0 * log(1 + exp(-alpha)) + (n1 + n2) * (log(1 + exp(alpha)));
  nllik = nllik + n1 * beta / 2 + sum1 * exp(-beta) / 2;
  nllik = nllik - sum(log(exp(y - C) + exp(alpha - C)) + C);

  Rcpp::NumericVector ret = Rcpp::as<Rcpp::NumericVector>(wrap(nllik));
  return ret;
}'

gradient.include <- 'Rcpp::NumericVector grad(SEXP xs, SEXP env) {
  Rcpp::NumericVector par(xs);
  Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);

  int fixpi0 = e["fixpi0"];
  int fixa = e["fixa"];
  int fixmu = e["fixmu"];

  double n0 = e["n0"];
  double n1 = e["n1"];
  double n2 = e["n2"];
  double sum1 = e["sum1"];
  arma::vec x = Rcpp::as<arma::vec>(e["x"]);
  arma::vec s2 = Rcpp::as<arma::vec>(e["s2"]);

  int i = 0;
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
  if (fixmu) {
    mu = e["mu"];
  } else {
    mu = par[i];
  }

  arma::vec z;
  if (fixmu) {
    z = Rcpp::as<arma::vec>(e["z"]);
  } else {
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
    grad[j] = -sum(tmp1 % (x - mu) / (1 + exp(alpha - y)) / s2);
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

likelihood.CPP <- inline::cxxfunction(signature(),
                                      body = likelihood.body,
                                      inc = likelihood.include,
                                      plugin = "RcppArmadillo")

gradient.CPP <- inline::cxxfunction(signature(),
                                    body = gradient.body,
                                    inc = gradient.include,
                                    plugin = "RcppArmadillo")
