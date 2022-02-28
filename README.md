# ebnm: Fit the empirical Bayes normal means problem

[![Travis Build Status](https://travis-ci.com/stephenslab/ebnm.svg?branch=master)](https://app.travis-ci.com/github/stephenslab/ebnm)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/l4u64gdn4noqlb1i?svg=true)](https://ci.appveyor.com/project/pcarbo/ebnm)
[![CircleCI build status](https://circleci.com/gh/stephenslab/ebnm.svg?style=svg)](https://app.circleci.com/pipelines/github/stephenslab/ebnm)
[![codecov](https://codecov.io/gh/stephenslab/ebnm/branch/master/graph/badge.svg)](https://app.codecov.io/gh/stephenslab/ebnm)

The `ebnm` package provides functions to solve the (heteroskedastic) "empirical Bayes normal means" (EBNM) problem for various choices of prior family. The model is

x_j | θ_j, s_j \sim N(θ_j, s_j^2)

θ_j | s_j \sim g \in G

where the distribution g is to be estimated. The distribution g is referred to as the "prior distribution" for θ and G is a specified family of prior distributions. Several options for G are implemented, some parametric and others non-parametric; see below for examples.


Solving the EBNM problem involves two steps. First, estimate g \in G via maximum marginal likelihood, yielding an estimate

\hat{g}:= \arg\max_{g \in G} L(g)

where

L(g):= ∏_j \int p(x_j | θ_j, s_j) g(dθ_j)

Second, compute the posterior distributions p(θ_j | x_j, s_j, \hat{g}) and/or summaries such as posterior means and posterior second moments.


The prior families that have been implemented include:

"point_normal":
The family of mixtures where one component is a point mass at μ and the other is a normal distribution centered at μ.

"point_laplace":
The family of mixtures where one component is a point mass at zero and the other is a double-exponential distribution.

"point_exponential":
The family of mixtures where one component is a point mass at zero and the other is a (nonnegative) exponential distribution.

"normal":
The family of normal distributions.

"horseshoe":
The family of horseshoe distributions.

"normal_scale_mixture":
The family of scale mixtures of normals.

"unimodal":
The family of all unimodal distributions.

"unimodal_symmetric":
The family of symmetric unimodal distributions.

"unimodal_nonnegative":
The family of unimodal distributions with support constrained to be greater than the mode.

"unimodal_nonpositive":
The family of unimodal distributions with support constrained to be less than the mode.

"npmle":
The family of all distributions.

"deconvolver":
A non-parametric exponential family with a natural spline basis. Like npmle, there is no unimodal assumption, but whereas npmle produces spiky estimates for g, deconvolver estimates are much more regular. See Narasimhan and Efron (2020) for details.

## License

The *ebnm* source code repository is free software: you can
redistribute it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html). All
the files in this project are part of *ebnm*. This project is
distributed in the hope that it will be useful, but **without any
warranty**; without even the implied warranty of **merchantability or
fitness for a particular purpose**.

## Quick Start

Install the ebnm package using `devtools`:

```R
library(devtools)
install_github("stephenslab/ebnm")
```

Load `ebnm` into your R environment, and get help:

```R
library(ebnm)
?ebnm
```

Try an example
```R
set.seed(1)
mu = c(rep(0,500),rnorm(500)) # true means
x = mu + rnorm(1000) # observations with standard error 1
x.ebnm = ebnm_point_normal(x, 1)
plot(mu, x.ebnm$posterior$mean) # plot posterior mean against true values
```

## Credits 

The `ebnm_point_laplace` function is derived from ideas and code in the *EbayesThresh* package
by I Johnstone and B Silverman.

