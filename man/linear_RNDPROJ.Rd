% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_RNDPROJ.R
\name{do.rndproj}
\alias{do.rndproj}
\title{Random Projection}
\usage{
do.rndproj(
  X,
  ndim = 2,
  preprocess = c("null", "center", "scale", "cscale", "whiten", "decorrelate"),
  type = c("gaussian", "achlioptas", "sparse"),
  s = max(sqrt(ncol(X)), 3)
)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations
and columns represent independent variables.}

\item{ndim}{an integer-valued target dimension.}

\item{preprocess}{an additional option for preprocessing the data.
Default is "null". See also \code{\link{aux.preprocess}} for more details.}

\item{type}{a type of random projection, one of "gaussian","achlioptas" or "sparse".}

\item{s}{a tuning parameter for determining values in projection matrix. While default
is to use \eqn{max(log \sqrt{p},3)}, it is required for \eqn{s \ge 3}.}
}
\value{
a named list containing
\describe{
\item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
\item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
\item{epsilon}{an estimated error \eqn{\epsilon} in accordance with JL lemma.}
\item{trfinfo}{a list containing information for out-of-sample prediction.}
}
}
\description{
\code{do.rndproj} is a linear dimensionality reduction method based on
random projection technique, featured by the celebrated Johnson–Lindenstrauss lemma.
}
\details{
The Johnson-Lindenstrauss(JL) lemma states that given \eqn{0 < \epsilon < 1}, for a set
\eqn{X} of \eqn{m} points in \eqn{R^N} and a number \eqn{n > 8log(m)/\epsilon^2},
there is a linear map \eqn{f:R^N} to {R^n} such that
\deqn{(1-\epsilon)|u-v|^2 \le |f(u)-f(v)|^2 \le (1+\epsilon)|u-v|^2}
for all \eqn{u,v} in \eqn{X}.

Three types of random projections are supported for an \code{(p-by-ndim)} projection matrix \eqn{R}.
\enumerate{
\item Conventional approach is to use normalized Gaussian random vectors sampled from unit sphere \eqn{S^{p-1}}.
\item Achlioptas suggested to employ a sparse approach using samples from \eqn{\sqrt{3}(1,0,-1)} with probability \eqn{(1/6,4/6,1/6)}.
\item Li et al proposed to  sample from \eqn{\sqrt{s}(1,0,-1)}
with probability \eqn{(1/2s,1-1/s,1/2s)} for \eqn{s\ge 3}
to incorporate sparsity while attaining speedup with little loss in accuracy. While
the original suggsetion from the authors is to use \eqn{\sqrt{p}} or \eqn{p/log(p)}
for \eqn{s}, any user-supported \eqn{s \ge 3} is allowed.
}
}
\examples{
## use iris data
data(iris)
set.seed(100)
subid = sample(1:150, 50)
X     = as.matrix(iris[subid,1:4])
label = as.factor(iris[subid,5])

## 1. Gaussian projection
output1 <- do.rndproj(X,ndim=2)

## 2. Achlioptas projection
output2 <- do.rndproj(X,ndim=2,type="achlioptas")

## 3. Sparse projection
output3 <- do.rndproj(X,type="sparse")

## Visualize three different projections
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(output1$Y, pch=19, col=label, main="RNDPROJ::Gaussian")
plot(output2$Y, pch=19, col=label, main="RNDPROJ::Arclioptas")
plot(output3$Y, pch=19, col=label, main="RNDPROJ::Sparse")
par(opar)

}
\references{
\insertRef{johnson_extensions_1984}{Rdimtools}

\insertRef{achlioptas_databasefriendly_2003}{Rdimtools}

\insertRef{li_very_2006}{Rdimtools}
}
\concept{linear_methods}
