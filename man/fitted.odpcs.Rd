\name{fitted.odpcs}
\alias{fitted.odpcs}
\title{Get Reconstructed Time Series From an odpcs Object}
\description{
  Get reconstructed time series from an odpcs object.
}


\usage{
  \method{fitted}{odpcs}(object, num_comp = 1, \dots)
}
\arguments{
  \item{object}{An object of class \code{odpcs}, usually the result of \link{odpc}.}
  \item{num_comp}{Integer indicating how many components to use for the reconstrucion. Default is 1.}
  \item{\dots}{Additional arguments for compatibility.}
}


\value{A matrix that is the reconstruction of the original series.}


\seealso{
  \code{\link{odpc}}, \code{\link{crit.odpc}}, \code{\link{cv.odpc}}
}

\examples{
  T <- 200 #length of series
  m <- 10 #number of series
  set.seed(1234)
  f <- rnorm(T + 1)
  x <- matrix(0, T, m)
  u <- matrix(rnorm(T * m), T, m)
  for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
  }
  fit <- odpc(x, ks = matrix(c(1, 1, 1, 0), 2, 2))
  recons <- fitted(fit, num_comp = 2)
}

\keyword{ts}