# ebnm: Fit the Empirical Bayes Normal Means (EBNM) problem

The `ebnm` package provides functions to solve the EBNM problem. 

Currently two functions are provided: `ebnm_point_normal` and `ebnm_point_laplace`
which solve the EBNM with a "point-normal" and "point-Laplace" prior respectively.
That is, the prior is a mixture
of a point mass at 0 and a normal or Laplace distribution. 

The goal is for the methods to be fast and numerically stable. 
The normal version is considerably faster (and probably more stable).


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
mu = c(rep(0,500),rnorm(500)) # true means
x = mu + rnorm(1000) #observations with standard error 1
x.ebnm = ebnm_point_normal(x,1)
plot(mu,ashr::get_pm(x.ebnm)) # plot posterior mean against true values
```

## Credits 

The `ebnm_point_laplace` function is derived from ideas and code in the *EbayesThresh* package
by I Johnstone and B Silverman.

