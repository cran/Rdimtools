% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/v1_feature_LSCORE.R
\name{do.lscore}
\alias{do.lscore}
\title{Laplacian Score}
\usage{
do.lscore(X, ndim = 2, ...)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations
and columns represent independent variables.}

\item{ndim}{an integer-valued target dimension (default: 2).}

\item{...}{extra parameters including \describe{
\item{preprocess}{an additional option for preprocessing the data.
See also \code{\link{aux.preprocess}} for more details (default: \code{"null"}).}
\item{type}{a vector of neighborhood graph construction. Following types are supported;
\code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
See also \code{\link{aux.graphnbd}} for more details (default: \code{c("proportion",0.1)}).}
\item{t}{bandwidth parameter for heat kernel in \eqn{(0,\infty)} (default: \code{1}).}
}}
}
\value{
a named \code{Rdimtools} S3 object containing
\describe{
\item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
\item{lscore}{a length-\eqn{p} vector of laplacian scores. Indices with smallest values are selected.}
\item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
\item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
\item{trfinfo}{a list containing information for out-of-sample prediction.}
\item{algorithm}{name of the algorithm.}
}
}
\description{
Laplacian Score \insertCite{he_laplacian_2005}{Rdimtools} is an unsupervised linear feature extraction method. For each
feature/variable, it computes Laplacian score based on an observation that data from the
same class are often close to each other. Its power of locality preserving property is used, and
the algorithm selects variables with smallest scores.
}
\examples{
\donttest{
## use iris data
## it is known that feature 3 and 4 are more important.
data(iris)
set.seed(100)
subid    <- sample(1:150, 50)
iris.dat <- as.matrix(iris[subid,1:4])
iris.lab <- as.factor(iris[subid,5])

## try different kernel bandwidth
out1 = do.lscore(iris.dat, t=0.1)
out2 = do.lscore(iris.dat, t=1)
out3 = do.lscore(iris.dat, t=10)

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(out1$Y, pch=19, col=iris.lab, main="bandwidth=0.1")
plot(out2$Y, pch=19, col=iris.lab, main="bandwidth=1")
plot(out3$Y, pch=19, col=iris.lab, main="bandwidth=10")
par(opar)
}

}
\references{
\insertAllCited{}
}
\author{
Kisung You
}
\concept{feature_methods}
