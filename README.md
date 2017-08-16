
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/esmucler/gdpc.svg?branch=master)](https://travis-ci.org/esmucler/gdpc)

odpc
====

This package provides functions for computing One-Sided Dynamic Principal Components, a novel multivariate time series dimension reduction technique proposed in Peña, Smucler and Yohai (2017).

------------------------------------------------------------------------

### Installation

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
f <- matrix(0, 2 * T, 1)
v <- rnorm(2 * T)
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
mse <- sum((x[T, ] - forecasts)**2)/(m)
mse
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
