# ebnm: Fit the empirical Bayes normal means problem

The `ebnm` package provides functions to solve the (heteroskedastic)
Empirical Bayes Normal Means (EBNM) problem for various choices of prior family.
The model is $$ x_j | \theta_j, s_j \sim N(\theta_j, s_j^2), $$
$$ \theta_j | s_j \sim g \in G, $$ where the distribution g (referred to as the 
"prior distribution" for $\theta$) is to be estimated and G is a specified family 
of prior distributions. Several options
for G are implemented, some parametric and others non-parametric.

Solving the EBNM problem involves
two steps. First, estimate g $\in$  G via maximum marginal likelihood,
yielding an estimate $$ \hat{g}:= \arg\max_{g \in G} L(g), $$ 
where $$ L(g):= \prod_j \int p(x_j | \theta_j, s_j)  g(d\theta_j). $$
Second, compute the posterior distributions 
$ p(\theta_j | x_j, s_j, \hat{g}) $ and/or summaries
such as posterior means and posterior second moments.

## License

The *ebnm* source code repository is free software: you can
redistribute it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html). All
the files in this project are part of *ebnm*. This project is
distributed in the hope that it will be useful, but **without any
warranty**; without even the implied warranty of **merchantability or
fitness for a particular purpose**. See file [LICENSE](LICENSE) for
the full text of the license.

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

