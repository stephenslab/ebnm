# ebnm: Fit the Empirical Bayes Normal Means problem

The `ebnm` package fits the EBNM problem. Currently it does this
for a simple "point-Laplace" prior: that is a mixture
of a point mass at 0 and a Laplace distribution. It is derived from the Empirical
Bayes methods developed by I. M. Johnstone and B. W. Silverman.
It is designed to be fast and stable.

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

There is only one exported function, `ebnm_point_laplace`.
Load `ebnm` into your R environment, and get help:

```R
library(ebnm)
??ebnm_point_laplace
```


## Credits 

The *ebnm* software package is based on code and idea from the *EbayesThresh* package.

