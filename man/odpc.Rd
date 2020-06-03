\name{odpc}
\alias{odpc}
\title{Fitting of One-Sided Dynamic Principal Components}
\description{
Computes One-Sided Dynamic Principal Components for a given number of lags.
}


\usage{odpc(Z, ks, method, tol = 1e-04, niter_max = 500)}
\arguments{
  \item{Z}{Data matrix. Each column is a different time series.}
  \item{ks}{Matrix or vector of integers. If a matrix, each row is the vector with number of lags to use for each              component. First column has the number of lags used to define the dynamic principal component (\eqn{k_1}), second column has the number of lags of the dynamic principal component used to reconstruct the series (\eqn{k_2}). If a vector, its entries are taken as both \eqn{k_1} and \eqn{k_2}             for each component}
  \item{method}{A string specifying the algorithm used. Options are 'ALS', 'mix' or 'gradient'. See details below.}
  \item{tol}{Relative precision. Default is 1e-4.}
  \item{niter_max}{Integer. Maximum number of iterations. Default is 500.}
}

\details{

Consider the vector time series \eqn{\mathbf{z}_{1},\dots,\mathbf{z}_{T}}, where
\eqn{\mathbf{z}_{t}=(z_{t,1},\dots,z_{t,m})^{\prime}}.  Let \eqn{\mathbf{a}%
=(\mathbf{a}_{0}^{\prime},\dots,\mathbf{a}_{k_{1}}^{\prime})^{\prime}}, where
\eqn{\mathbf{a}_{h}^{\prime}=(a_{h,1},...,a_{h,m})}, be a vector of dimension
\eqn{m(k_{1}+1)\times1}, let \eqn{\boldsymbol{\alpha}^{\prime}=(\alpha_{1}%
,\dots,\alpha_{m})} and \eqn{\mathbf{B}} the matrix that has coefficients
\eqn{b_{h,j}} and dimension \eqn{(k_{2}+1)\times m}. Consider
\deqn{
f_{t}=\sum\limits_{j=1}^{m}\sum\limits_{h=0}^{k_{1}}a_{h,j}z_{t-h,j}\quad
t=k_{1}+1,\dots,T, \nonumber
}
and suppose we use \eqn{f_t} and \eqn{k_{2}} of its lags to reconstruct the series
as
\deqn{
  z_{t,j}^{R}(\mathbf{a},\boldsymbol{\alpha},\mathbf{B)}=\alpha_{j}
               +\sum\limits_{h=0}^{k_{2}}b_{h,j}f_{t-h}.\nonumber
}
Let
\deqn{
MSE(\mathbf{a},\boldsymbol{\alpha},\mathbf{B})=\frac{1}{T-(k_{1}%
+k_{2})}\sum\limits_{j=1}^{m}\sum\limits_{t=(k_{1}+k_{2}%
)+1}^{T}(z_{t,j}-z_{t,j}^{R}(\mathbf{a},\boldsymbol{\alpha},\mathbf{B)})^{2}
}
be the reconstruction MSE.
The first one-sided dynamic principal component is defined as the series
\deqn{
\widehat{f}_{t}=\sum\limits_{j=1}^{m}\sum\limits_{h=0}^{k_{1}}\widehat{a}_{h,j}z_{t-h,j}\quad
t=k_{1}+1,\dots,T, \nonumber
}
for optimal values \eqn{(\widehat{\mathbf{a}},\widehat{\boldsymbol{\alpha}}%
              ,\widehat{\mathbf{B}})}
that satisfy
\deqn{
  MSE(\widehat{\mathbf{a}},\widehat{\boldsymbol{\alpha}}%
              ,\widehat{\mathbf{B}})=\min_{\Vert\mathbf{a}\Vert=1,\boldsymbol{\alpha
              },\mathbf{B}}MSE(\mathbf{a},\boldsymbol{\alpha},\mathbf{B}).
  \nonumber
}

The second one-sided dynamic principal component is defined similarly, but now the
residuals of the first one-sided dynamic principal component are to be reconstructed.

If method = 'ALS', an Alternating Least Squares type algorithm is used to compute the solution. If 'mix' is chosen, in each iteration Least Squares is used to compute the matrix of loadings and intercepts, but one iteration of Coordinate Descent is performed to compute the vector a that defines the dynamic principal component. If method = 'gradient', in each iteration Least Squares is used to compute the matrix of loadings and intercepts, but one iteration of Gradient Descent is performed to compute the vector a that defines the dynamic principal component. By default, 'ALS' is used when the number of series is less than 40, else 'mix' is used.
}

\value{
  An object of class \code{odpcs}, that is, a list of length equal to the number of computed components. The i-th entry of this list is an object of class \code{odpc}, that is, a list with entries
  \item{f}{Coordinates of the i-th dynamic principal component corresponding to the periods \eqn{k_1 + 1,\dots,T}.}
  \item{mse}{Mean squared error of the reconstruction using the first i components.}
  \item{k1}{Number of lags used to define the i-th dynamic principal component f.}
  \item{k2}{Number of lags of f used to reconstruct.}
  \item{alpha}{Vector of intercepts corresponding to f.}
  \item{a}{Vector that defines the i-th dynamic principal component}
  \item{B}{Matrix of loadings corresponding to f. Row number \eqn{k} is the vector of \eqn{k-1} lag loadings.}
      \item{call}{The matched call.}
      \item{conv}{Logical. Did the iterations converge?}
  \code{components}, \code{fitted}, \code{plot} and \code{print} methods are available for this class.
}

\references{
Peña D., Smucler E. and Yohai V.J. (2019). “Forecasting Multiple Time Series with One-Sided Dynamic Principal Components.” Journal of the American Statistical Association.
}

\seealso{
\code{\link{crit.odpc}}, \code{\link{cv.odpc}}, \code{\link{plot.odpc}}, \code{\link{fitted.odpcs}}, \code{\link{components_odpcs}}, \code{\link{forecast.odpcs}}
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
fit
}
\keyword{ts}