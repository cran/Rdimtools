#' Neighborhood Preserving Embedding
#'
#' \code{do.npe} performs a linear dimensionality reduction using Neighborhood Preserving
#' Embedding (NPE) proposed by He et al (2005). It can be regarded as a linear approximation
#' to Laplacian Eigenmaps.
#'
#' Like Locally Linear Embedding (LLE), it is possible for the weight matrix being rank deficient.
#' If \code{regtype} is set to \code{TRUE} with a proper value of \code{regparam}, it will
#' perform Tikhonov regularization as designated. When regularization is needed
#' with \code{regtype} parameter to be \code{FALSE}, it will automatically find a suitable
#' regularization parameter and put penalty for stable computation. See also
#' \code{\link{do.lle}} for more details.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param weight \code{TRUE} to perform NPE on weighted graph, or \code{FALSE} otherwise.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null" and three options of "center", "decorrelate", or "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param regtype \code{FALSE} for not applying automatic Tikhonov Regularization,
#' or \code{TRUE} otherwise.
#' @param regparam a positive real number for Regularization. Default value is 1.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \code{(n-by-ndim)} matrix whose rows are embedded observations.}
#' \item{eigval}{a vector of eigenvalues corresponding to basis expansion in an ascending order.}
#' \item{projection}{a \code{(p-by-ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'@examples
#'# generate data
#'X <- aux.gensamples(n=333)
#'
#'## 1. connecting 5% of data for graph construction.
#'output1 <- do.npe(X,ndim=2,type=c("proportion",0.05))
#'
#'## 2. constructing 25%-connected graph with regularization parameter
#'output2 <- do.npe(X,ndim=2,type=c("proportion",0.2),symmetric='intersect',regparam=1.0)
#'
#'## 3. constructing half-connected graph with reg parameter = 10.0.
#'output3 <- do.npe(X,ndim=2,type=c("proportion",0.5),regparam=10.0)
#'
#'## Visualize three different projections
#'par(mfrow=c(1,3))
#'plot(output1$Y[,1],output1$Y[,2],main="5%")
#'plot(output2$Y[,1],output2$Y[,2],main="25%")
#'plot(output3$Y[,1],output3$Y[,2],main="50%")
#'
#' @references
#' \insertRef{he_neighborhood_2005}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_NPE
#' @export
do.npe <- function(X,ndim=2,type=c("proportion",0.1),symmetric="union",weight=TRUE,preprocess="null",regtype=FALSE,regparam=1){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.npe : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. aux.graphnbd
  #   type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   symmetric : 'intersect','union', or 'asymmetric'
  # 2-2. NPE itself
  #   weight     : TRUE
  #   preprocess : 'null', 'center','decorrelate', or 'whiten'
  #   regtype    : FALSE (default; need a value) / TRUE
  #   regparam   : 1 (default)

  nbdtype = type
  nbdsymmetric = symmetric
  algweight = weight
  algpreprocess = preprocess
  if (!is.element(algpreprocess,c("null","center","whiten","decorrelate"))){
    stop("* do.npe : 'preprocess' argument is invalid.")
  }
  if (!is.logical(regtype)){stop("* do.npe : 'regtype' should be a logical variable.")}
  if (!is.numeric(regparam)||is.na(regparam)||is.infinite(regparam)||(regparam<=0)){
    stop("* do.npe : 'regparam' should be a positive real-valued number; it is a Tikhonov Regularization Factor.")}


  #   regtype    : FALSE (default; need a value) / TRUE
  #   regparam   : 1 (default)

  # 3. process : data preprocessing
  if (algpreprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X))
  } else {
    tmplist = aux.preprocess(X,type=algpreprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }

  n = nrow(pX)
  p = ncol(pX)

  # 4. process : neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)


  # 5. main 1 : compute Weights
  #   k = max(apply(nbdstruct$mask,2,function(x) sum(as.double((x==TRUE)))))
  W = array(0,c(n,n))
  if (regtype==TRUE){
    regvals = array(0,c(1,n))
  }
  for (i in 1:n){
    #   5-1. separate target mask vector
    tgtmask = nbdstruct$mask[i,]
    tgtidx  = which(tgtmask==TRUE)

    #   5-2. select data
    #        For convenience, target matrix is transposed for Armadillo
    vec_tgt = pX[i,]
    mat_tgt = t(pX[tgtidx,])
    k = ncol(mat_tgt)

    #   5-3. compute with regularization
    #   5-3-1. No Automatic Regularization
    if (regtype==FALSE){
      w = method_lleW(mat_tgt,vec_tgt,regparam);
    } else {
      #   5-3-2. Automatic Regularization
      outW = method_lleWauto(mat_tgt,vec_tgt);
      w          = outW$w
      regvals[i] = outW$regparam
    }
    W[i,tgtidx] = w;
  }

  # 6. main computation
  tpX = t(pX)
  output = method_npe(tpX,W);

  eigvals = output$eigval
  eigvecs = output$eigvec

  # 7. return output
  result = list()
  result$Y = pX %*% eigvecs[,1:ndim]
  result$eigval = eigvals
  result$projection = eigvecs[,1:ndim]
  trfinfo$algtype   = "linear"
  result$trfinfo = trfinfo
  return(result)
}