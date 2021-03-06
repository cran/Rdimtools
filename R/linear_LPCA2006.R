#' Locally Principal Component Analysis by Yang et al. (2006)
#'
#' Locally Principal Component Analysis (LPCA) is an unsupervised linear dimension reduction method.
#' It focuses on the information brought by local neighborhood structure and seeks the corresponding
#' structure, which may contain useful information for revealing discriminative information of the data.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris dataset
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,100)
#' X     = as.matrix(iris[subid,1:4])
#' lab   = as.factor(iris[subid,5])
#'
#' ## try different neighborhood size
#' out1 <- do.lpca2006(X, ndim=2, type=c("proportion",0.25))
#' out2 <- do.lpca2006(X, ndim=2, type=c("proportion",0.50))
#' out3 <- do.lpca2006(X, ndim=2, type=c("proportion",0.75))
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=lab, main="LPCA2006::25% connected")
#' plot(out2$Y, pch=19, col=lab, main="LPCA2006::50% connected")
#' plot(out3$Y, pch=19, col=lab, main="LPCA2006::75% connected")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{yang_locally_2006}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LPCA2006
#' @concept linear_methods
#' @export
do.lpca2006 <- function(X, ndim=2, type=c("proportion",0.1),
                    preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lpca2006 : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. neighborhood information : asymmetric
  nbdtype = type
  nbdsymmetric = "asymmetric"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR LOCALLY PCA
  #   1.  build H and symmetrize
  H = nbdmask*1.0
  H = (H + t(H))/2
  #   2. L for graph laplacian
  L = base::diag(base::rowSums(H))-H
  #   3. spectral decomposition of L
  Pl = lpca_spectralhalf(L)
  #   4. compute R
  R  = t(Pl)%*%pX%*%t(pX)%*%Pl
  #   5. compute top eigenpairs for R
  topR = RSpectra::eigs(R, ndim)
  #   6. projection matrix as noted with scaling
  projection = aux.adjprojection((t(pX)%*%Pl%*%(topR$vectors))%*%diag(1/sqrt(topR$values)))

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}


#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
lpca_spectralhalf <- function(A){
  eigA = base::eigen(A)

  Aval = sqrt(eigA$values)
  Pl   = (eigA$vectors)%*%diag(Aval)
  return(Pl)
}
