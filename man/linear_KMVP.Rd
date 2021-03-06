% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_KMVP.R
\name{do.kmvp}
\alias{do.kmvp}
\title{Kernel-Weighted Maximum Variance Projection}
\usage{
do.kmvp(
  X,
  label,
  ndim = 2,
  preprocess = c("center", "scale", "cscale", "decorrelate", "whiten"),
  bandwidth = 1
)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations
and columns represent independent variables.}

\item{label}{a length-\eqn{n} vector of data class labels.}

\item{ndim}{an integer-valued target dimension.}

\item{preprocess}{an additional option for preprocessing the data.
Default is "center". See also \code{\link{aux.preprocess}} for more details.}

\item{bandwidth}{bandwidth parameter for heat kernel as the equation above.}
}
\value{
a named list containing
\describe{
\item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
\item{trfinfo}{a list containing information for out-of-sample prediction.}
\item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
}
}
\description{
Kernel-Weighted Maximum Variance Projection (KMVP) is a generalization of
Maximum Variance Projection (MVP). Even though its name contains \emph{kernel}, it is
not related to kernel trick well known in the machine learning community. Rather, it
generalizes the binary penalization on class discrepancy,
\deqn{S_{ij} = \exp(-\|x_i-x_j\|^2/t) \quad\textrm{if}\quad C_i \ne C_j}
where \eqn{x_i} is an \eqn{i}-th data point and \eqn{t} a kernel bandwidth (\code{bandwidth}). \bold{Note} that
when the bandwidth value is too small, it might suffer from numerical instability and rank deficiency due to its formulation.
}
\examples{
## use iris data
data(iris)
set.seed(100)
subid = sample(1:150, 50)
X     = as.matrix(iris[subid,1:4])
label = as.factor(iris[subid,5])

## perform KMVP with different bandwidths
out1 = do.kmvp(X, label, bandwidth=0.1)
out2 = do.kmvp(X, label, bandwidth=1)
out3 = do.kmvp(X, label, bandwidth=10)

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(out1$Y, main="bandwidth=0.1", col=label, pch=19)
plot(out2$Y, main="bandwidth=1",   col=label, pch=19)
plot(out3$Y, main="bandwidth=10",  col=label, pch=19)
par(opar)

}
\references{
\insertRef{zhang_maximum_2007}{Rdimtools}
}
\seealso{
\code{\link{do.mvp}}
}
\author{
Kisung You
}
\concept{linear_methods}
