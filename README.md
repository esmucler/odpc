
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/esmucler/odpc.svg?branch=master)](https://travis-ci.org/esmucler/odpc) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/odpc)](https://cran.r-project.org/package=odpc) [![Downloads](http://cranlogs.r-pkg.org/badges/odpc)](https://cran.r-project.org/package=odpc)

odpc
====

This package provides functions for computing One-Sided Dynamic Principal Components, a novel multivariate time series dimension reduction technique proposed in [Peña, Smucler and Yohai (2019)](https://doi.org/10.1080/01621459.2018.1520117).

Version 2.0.4 includes a new feature, not discussed in the original paper, that of fitting sparse one-sided dynamic principal components.
These are constructed by penalizing the L1 norm of the vector of coefficients that defines the one-sided dynamic principal components. 
This new loss function is minimized using a new algorithm based on alternating proximal gradient descent. The penalty parameters are chosen automatically by minimizing a BIC type criterion. See the crit.sparse_odpc function. The original one-sided dynamic principal components (i.e., the non-sparse ones) are now also optimized using alternating gradient descent by default. Earlier algorithms based on alternating least squares and coordinate descent optimization are still available as possible options for the method argument. 

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=odpc).

``` r
install.packages('odpc', dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/esmucler/odpc)

``` r
library(devtools)
devtools::install_github("esmucler/odpc")
```

### Usage

``` r
library(odpc)
# An example using an artificial data set
T <- 201 #length of series
m <- 10 #number of series
set.seed(1234)
f <- matrix(0, 2 * T + 1, 1)
v <- rnorm(2 * T + 1)
f[1] <- rnorm(1)
theta <- 0.7
# f satisfies and AR(1) model
for (t in  2:(2 * T)){
  f[t] <- theta * f[t - 1] + v[t]
}
f <- f[T:(2 * T)]
x <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
  x[, i] <- sin(2 * pi * (i/m)) * f[1:T] + cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
}
# Compute one odpc with one lag
fit <- odpc(x[1:(T - 1), ], ks = c(1))
# Get a one step ahed forecast of x
forecasts <- forecast.odpcs(fit, h = 1)
mse <- mean((x[T, ] - forecasts)**2)
mse
# Compute one sparse odpc with one lag
sparse_fit <- crit.sparse_odpc(x[1:(T - 1), ], k_list = c(1))
sparse_fit[[1]]$a
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
