\name{components.odpcs}
\alias{components.odpcs}
\title{Get One-Sided Dynamic Principal Components From an odpcs Object}
\description{
Get One-Sided Dynamic Principal Components from an \code{odpcs} object.
}


\usage{
\method{components}{odpcs}(object, which_comp = 1)
}
\arguments{
  \item{object}{An object of class \code{odpcs}, usually the result of \link{odpc}.}
  \item{which_comp}{Numeric vector indicating which components to get. Default is 1.}
}


\value{A list whose entries are the desired dynamic principal components.}
\author{
Daniel Peña, Ezequiel Smucler, Victor Yohai
}



\seealso{
\code{\link{odpc}}, \code{\link{forecast.odpcs}}
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
  comps <- components.odpcs(fit, which_comp = c(1, 2))
}

\keyword{ts}