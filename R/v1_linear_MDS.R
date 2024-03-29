#' (Classical) Multidimensional Scaling
#'
#' \code{do.mds} performs a classical Multidimensional Scaling (MDS) using
#' \code{Rcpp} and \code{RcppArmadillo} package to achieve faster performance than
#' \code{\link[stats]{cmdscale}}.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param ... extra parameters including \describe{
#' \item{preprocess}{an additional option for preprocessing the data.
#' Default is \code{"null"}. See also \code{\link{aux.preprocess}} for more details.}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' lab   = as.factor(iris[subid,5])
#'
#' ## compare with PCA
#' Rmds <- do.mds(X, ndim=2)
#' Rpca <- do.pca(X, ndim=2)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(Rmds$Y, pch=19, col=lab, main="MDS")
#' plot(Rpca$Y, pch=19, col=lab, main="PCA")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{kruskal_multidimensional_1964}{Rdimtools}
#'
#' @concept linear_methods
#' @rdname linear_MDS
#' @export
do.mds <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  ## BASIC
  #   explicit
  if (!is.matrix(X)){stop("* do.mds : 'X' should be a matrix.")}
  myndim = min(max(1, round(ndim)), ncol(X)-1)

  #   implicit
  params = list(...)
  pnames = names(params)
  if ("preprocess"%in%pnames){
    par_preprocess = tolower(params$preprocess)
  } else {
    par_preprocess = "null"
  }

  #------------------------------------------------------------------------
  ## COMPUTE
  #   preprocessing
  tmplist = aux.preprocess.hidden(X, type=par_preprocess, algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   compute
  output = dt_mds(pX, myndim)
  output$trfinfo = trfinfo
  return(structure(output, class="Rdimtools"))
}

# call for later use ------------------------------------------------------
#' @keywords internal
pydo_mds <- function(myX, mydim){
  return(do.mds(myX, ndim=mydim))
}
