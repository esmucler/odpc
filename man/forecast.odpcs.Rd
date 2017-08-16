\name{forecast.odpcs}
\alias{forecast.odpcs}
\title{Get Forecast From an odpcs Object}
\description{
  Get forecasts from an \code{odpcs} object.
}


\usage{
  \method{forecast}{odpcs}(object, h, Z = NULL, add_residuals = FALSE, \dots)
}
\arguments{
  \item{object}{An object of class \code{odpcs}, usually the result of \link{odpc}.}
  \item{h}{Integer. Number of periods for forecasting.}
  \item{Z}{Original data. Only used if add_residuals = TRUE.}
  \item{add_residuals}{Logical. Should the forecasts of the reconstruction residuals
                       be added to the final forecast? Default is FALSE.}
  \item{\dots}{Additional arguments to be passed to \code{\link[forecast]{auto.arima}}.}
}

\details{

Suppose \eqn{q} dynamic principal components were fitted to the data, each with
\eqn{(k_{1}^{i},k_{2}^{i})} lags, \eqn{i=1,\dots,q}. Let \eqn{\widehat{\mathbf{f}}_{T}%
^{i}} be the vector with the estimated values for the \eqn{i}-th dynamic principal
component and \eqn{\widehat{\mathbf{B}}^{i}}, \eqn{\widehat{\boldsymbol{\alpha}}^{i}}
be the corresponding loadings and intercepts. 
Forecasts of the series are built by first fitting a SARIMA model to the components
using \code{\link[forecast]{auto.arima}} and getting their forecasts using \code{\link[forecast]{forecast.Arima}}.
Let \eqn{\widehat{f}_{T+h|T}^{i}} for \eqn{h>0} be the forecast of \eqn{f_{T+h}^{i}}
with information until time \eqn{T}. Then the \eqn{h}-steps ahead forecast of
\eqn{\mathbf{z}_{T}} is obtained as
\deqn{
\widehat{z}_{T+h|T,j}=\sum\limits_{i=1}^{q}\left(  \widehat{\alpha}_{j}%
^{i}+\sum\limits_{v=0}^{k_{2}^{i}}\widehat{b}_{v,j}^{i}\widehat{f}%
_{T+h-v|T}^{i}\right)  \quad j=1,\dots,m.
}

If add_residuals = TRUE, univariate SARIMA models are fitted to the residuals of the reconstruction, and their
forecasts are added to the forecasts described above.
}


\value{ A matrix that is the h-step ahead forecast of the original series.}
\author{
  Daniel Peña, Ezequiel Smucler, Victor Yohai
}



\seealso{
  \code{\link{odpc}}, \code{\link{components.odpcs}}, \code{\link{fitted.odpcs}}, \code{\link[forecast]{auto.arima}}, \code{\link[forecast]{forecast.Arima}}
}

\examples{
  T <- 201 #length of series
  m <- 10 #number of series
  set.seed(1234)
  f <- matrix(0, 2 * T, 1)
  v <- rnorm(2 * T)
  f[1] <- rnorm(1)
  theta <- 0.7
  for (t in  2:(2 * T)){
    f[t] <- theta * f[t - 1] + v[t]
  }
  f <- f[T:(2 * T)]
  x <- matrix(0, T, m)
  u <- matrix(rnorm(T * m), T, m)
  for (i in 1:m) {
    x[, i] <- sin(2 * pi * (i/m)) * f[1:T] + cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
  }
  fit <- odpc(x[1:(T - 1), ], ks = c(1))
  forecasts <- forecast.odpcs(fit, h = 1)
  mse <- mean((x[T, ] - forecasts)**2)
  mse
}

\keyword{ts}