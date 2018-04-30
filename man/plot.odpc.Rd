\name{plot.odpc}
\alias{plot.odpc}
\title{Plot One-Sided Dynamic Principal Components}
\description{
Plots an \code{odpc} object.
}


\usage{
\method{plot}{odpc}(x, which = 'Component', which_load = 0, \dots)
}
\arguments{
  \item{x}{An object of class \code{odpc}, usually one of the entries of the result of \link{odpc}.}
  \item{which}{String. Indicates what to plot, either 'Component' or 'Loadings'. Default is 'Component'.}
  \item{which_load}{Lag number indicating which loadings should be plotted. Only used if which = 'Loadings'. Default is 0.}
  \item{\dots}{Additional arguments to be passed to the plotting functions.}
}



\seealso{
\code{\link{odpc}}
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
fit <- odpc(x, ks = c(1))
plot(fit[[1]], xlab = '', ylab = '')
}

\keyword{ts}