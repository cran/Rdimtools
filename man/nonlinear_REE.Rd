% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonlinear_REE.R
\name{do.ree}
\alias{do.ree}
\title{Robust Euclidean Embedding}
\usage{
do.ree(
  X,
  ndim = 2,
  W = NA,
  preprocess = c("null", "center", "scale", "cscale", "whiten", "decorrelate"),
  initc = 1,
  dmethod = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
  maxiter = 100,
  abstol = 0.001
)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations
and columns represent independent variables.}

\item{ndim}{an integer-valued target dimension.}

\item{W}{an \eqn{(n\times n)} weight matrix. Default is uniform weight of 1s.}

\item{preprocess}{an additional option for preprocessing the data.
Default is "null". See also \code{\link{aux.preprocess}} for more details.}

\item{initc}{initial \code{c} value for subgradient iterating stepsize, \eqn{c/\sqrt{i}}.}

\item{dmethod}{a type of distance measure. See \code{\link[stats]{dist}} for more details.}

\item{maxiter}{maximum number of iterations for subgradient descent method.}

\item{abstol}{stopping criterion for subgradient descent method.}
}
\value{
a named list containing
\describe{
\item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
\item{niter}{the number of iterations taken til convergence. }
\item{trfinfo}{a list containing information for out-of-sample prediction.}
}
}
\description{
Robust Euclidean Embedding (REE) is an embedding procedure exploiting
robustness of \eqn{\ell_1} cost function. In our implementation, we adopted
a generalized version with weight matrix to be applied as well. Its original
paper introduced a subgradient algorithm to overcome memory-intensive nature of
original semidefinite programming formulation.
}
\examples{
\donttest{
## use iris data
data(iris)
set.seed(100)
subid = sample(1:150,50)
X     = as.matrix(iris[subid,1:4])
label = as.factor(iris[subid,5])

## try different distance method
output1 <- do.ree(X, maxiter=50, dmethod="euclidean")
output2 <- do.ree(X, maxiter=50, dmethod="maximum")
output3 <- do.ree(X, maxiter=50, dmethod="canberra")

## visualize three different projections
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(output1$Y, col=label, pch=19, main="dmethod-euclidean")
plot(output2$Y, col=label, pch=19, main="dmethod-maximum")
plot(output3$Y, col=label, pch=19, main="dmethod-canberra")
par(opar)
}

}
\references{
\insertRef{cayton_robust_2006}{Rdimtools}
}
\author{
Kisung You
}
\concept{nonlinear_methods}
