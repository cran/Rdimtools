% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonlinear_NNP.R
\name{do.nnp}
\alias{do.nnp}
\title{Nearest Neighbor Projection}
\usage{
do.nnp(
  X,
  ndim = 2,
  preprocess = c("null", "center", "scale", "cscale", "whiten", "decorrelate")
)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations
and columns represent independent variables.}

\item{ndim}{an integer-valued target dimension.}

\item{preprocess}{an additional option for preprocessing the data.
Default is "null". See also \code{\link{aux.preprocess}} for more details.}
}
\value{
a named list containing
\describe{
\item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
\item{trfinfo}{a list containing information for out-of-sample prediction.}
}
}
\description{
Nearest Neighbor Projection is an iterative method for visualizing high-dimensional dataset
in that a data is sequentially located in the low-dimensional space by maintaining
the triangular distance spread of target data with its two nearest neighbors in the high-dimensional space.
We extended the original method to be applied for arbitrarily low-dimensional space. Due the generalization,
we opted for a global optimization method of \emph{Differential Evolution} (\code{\link[RcppDE]{DEoptim}}) within in that it may add computational burden to certain degrees.
}
\examples{
\donttest{
## use iris data
data(iris)
set.seed(100)
subid = sample(1:150,50)
X     = as.matrix(iris[subid,1:4])
label = as.factor(iris[subid,5])

## let's compare with other methods
out1 <- do.nnp(X, ndim=2)      # NNP
out2 <- do.pca(X, ndim=2)      # PCA
out3 <- do.dm(X, ndim=2)     # Diffusion Maps

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(out1$Y, pch=19, col=label, main="NNP")
plot(out2$Y, pch=19, col=label, main="PCA")
plot(out3$Y, pch=19, col=label, main="Diffusion Maps")
par(opar)
}

}
\references{
\insertRef{tejada_improved_2003}{Rdimtools}
}
\author{
Kisung You
}
\concept{nonlinear_methods}
