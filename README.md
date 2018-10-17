# ebnm: Fit the Empirical Bayes Normal Means (EBNM) problem

The `ebnm` package provides functions to solve the (heteroskedastic)
Empirical Bayes Normal Means (EBNM) problem, which is as follows.
Observations $x=(x_1,\dots,x_n)$ are assumed to be independent with
$$x_j | \theta_j \sim N(\theta_j, s_j^2)$$ where the standard
deviations $s_j$ are assumed known and the means $\theta_j$ are to be
estimated.

In addition, the $\theta_j$ are assumed to be independent and identically distributed,
$$\theta_j \sim g() \in G$$ for some pre-specified family $G$.

Solving the EBNM problem involves two steps:

	- estimate $g$ by maximizing the likelihood. That is, $\hat{g} = \arg\max_{g \in G} p(x_1,\dots,x_n | g)$.

	- compute (summaries of) the posterior distributions $p(\theta_j | \hat{g}, x)$.
	
Currently two functions are provided: `ebnm_point_normal` and `ebnm_point_laplace`
which solve the EBNM with a "point-normal" and "point-Laplace" prior respectively.

That is, `ebnm_point_normal` solves EBNM with 
$$g(\cdot) = \pi_0 \delta_0(\cdot) + (1-\pi_0) N(\cdot; 0,1/a)$$ 
where $\delta_0$ denotes a point mass at 0, and $N(\cdot; \mu,\sigma^2)$ denotes
the density of a normal distribution with mean $\mu$ and variance $\sigma^2$. (So $a$ is the inverse-variance, or precision.)

And `ebnm_point_laplacel` solves EBNM with 
$$g(\cdot) = \pi_0 \delta_0(\cdot) + (1-\pi_0) DExp(\cdot; a)$$ 
where $\delta_0$ denotes a point mass at 0, and $DExp(\cdot; \lambda)$ denotes
the density of a double exponential (Laplace) distribution with rate parameter $\lamdba$.

In both cases the parameters $\pi_0$ and $a$ are to be estimated.

The goal is for the methods to be fast and numerically stable. 
The `ebnm_point_normal` function is considerably faster (and probably more stable).


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
??ebnm_point_normal
??ebnm_point_laplace
```

Try an example
```R
set.seed(1)
mu = c(rep(0,500),rnorm(500)) # true means
x = mu + rnorm(1000) #observations with standard error 1
x.ebnm = ebnm_point_normal(x,1)
plot(mu,ashr::get_pm(x.ebnm)) # plot posterior mean against true values
```

## How to build the webpages

Run the following commands in R from the [vignettes](vignettes)
directory:

```R
library(rmarkdown)
render("opt_normal.Rmd",output_dir = "../docs")
```

## Credits 

The `ebnm_point_laplace` function is derived from ideas and code in the *EbayesThresh* package
by I Johnstone and B Silverman.

