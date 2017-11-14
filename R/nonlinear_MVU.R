#' Maximum Variance Unfolding / Semidefinite Embedding
#'
#' The method of Maximum Variance Unfolding(MVU), also known as Semidefinite Embedding(SDE) is, as its names suggest,
#' to exploit semidefinite programming in performing nonlinear dimensionality reduction by \emph{unfolding}
#' neighborhood graph constructed in the original high-dimensional space. Its unfolding generates a gram
#' matrix \eqn{K} in that we can choose from either directly finding embeddings (\code{"spectral"}) or
#' use again Kernel PCA technique (\code{"kpca"}) to find low-dimensional representations. Note that
#' since \code{do.mvu} depends on \url{https://CRAN.R-project.org/package=Rcsdp}{Rcsdp}, we cannot guarantee its computational
#' efficiency as we are given a large dataset.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is ``null'', and other methods of ``decorrelate'',``center'' , and ``whiten'' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param projtype type of method for projection; either \code{"spectral"} or \code{"kpca"} used.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \code{(n-by-ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## generate ribbon-shaped data
#' ## in order to pass CRAN pretest, n is set to be small.
#' X = aux.gensamples(dname="ribbon",n=28)
#'
#' ## 1. standard MVU
#' output1 <- do.mvu(X,ndim=2)
#'
#' ## 2. standard setting with "kpca"-type projection
#' output2 <- do.mvu(X,ndim=2,projtype="kpca")
#'
#' ## 3. standard MVU for densly connected graph
#' output3 <- do.mvu(X,ndim=2,type=c("proportion",0.5))
#'
#' ## Visualize three different projections
#' par(mfrow=c(1,3))
#' plot(output1$Y[,1],output1$Y[,2],main="standard")
#' plot(output2$Y[,1],output2$Y[,2],main="kpca projection")
#' plot(output3$Y[,1],output3$Y[,2],main="densely connected graph")
#'
#' @references
#' \insertRef{weinberger_unsupervised_2006}{Rdimtools}
#'
#' @author Kisung You
#' @aliases do.sde
#' @rdname nonlinear_MVU
#' @export
do.mvu <- function(X,ndim=2,type=c("proportion",0.1),preprocess="null",projtype="spectral"){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.mvu : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. aux.graphnbd
  #   type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  # 2-2. mvu itself
  #   preprocess : 'center','decorrelate', 'null', or 'whiten'
  nbdtype = type
  nbdsymmetric = "union"

  algpreprocess = preprocess
  if (!is.element(algpreprocess,c("null","center","whiten","decorrelate"))){
    stop("* do.mvu : 'preprocess' argument is invalid.")
  }

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
  trfinfo$algtype = "nonlinear"

  n = nrow(pX)
  p = ncol(pX)

  # 4. process : neighborhood selection + ALPHA for SDE/MVU
  #            : make a triangle.
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric="union")
  M1 = nbdstruct$mask
  M2 = M1
  for (it in 1:nrow(M1)){
    for (i in 1:ncol(M1)){
      for (j in 1:ncol(M1)){
        if ((M1[it,i]==TRUE)&&(M1[it,j]==TRUE)){
          M2[i,j] = TRUE
          M2[j,i] = TRUE
        }
      }
    }
  }
  diag(M2) = FALSE

  # 5. Compute K
  #   5-1. basic settings
  N  = ncol(M2) # number of data points
  nconstraints = sum(M2)/2
  D2 = ((as.matrix(dist(pX)))^2)
  C  = list(.simple_triplet_diag_sym_matrix(1,N))
  K  = list(type="s",size=N)
  #   5-2. iterative computation
  iter = 1
  A = list()
  b = c()
  for (it1 in 1:nrow(M2)){
    for (it2 in it1:ncol(M2)){
      if (M2[it1,it2]==TRUE){
        tmpA = list(simple_triplet_sym_matrix(i=c(it1,it1,it2,it2),j=c(it1,it2,it1,it2),v=c(1,-1,-1,1),n=N))
        tmpb = D2[it1,it2]
        A[[iter]] = tmpA
        b[iter]   = tmpb
        iter = iter+1
      }
    }
  }
  #   5-3. final update for mean centered constraint
  A[[iter]] = list(matrix(1,N,N)/N)
  b[iter]   = 0
  outCSDP = csdp(C,A,b,K)
  KK      = as.matrix(outCSDP$X[[1]])

  projtype
  if (projtype=="spectral"){
    # 6. Embedding : Spectral Method, directly from K
    KKeigen = eigen(KK)
    eigvals = KKeigen$values
    eigvecs = KKeigen$vectors

    Y = (diag(sqrt(eigvals[1:ndim])) %*% t(eigvecs[,1:ndim]))
  } else if (projtype=="kpca"){
    # 7. Embedding : Kernel PCA method
    onesN = matrix(1,n,n)
    KKcentered = KK - (onesN%*%KK/n) - (KK%*%onesN/n) + (onesN%*%KK%*%onesN/(n^2))
    KKceigen   = eigen(KKcentered)

    Y     = (t(KKceigen$vectors[,1:ndim]) %*% KK)
  } else {
    stop("* do.mvu : 'projtype' should be either 'spectral' or 'kpca'.")
  }

  # 8. Return output
  result = list()
  result$Y = t(Y)
  result$trfinfo = trfinfo
  return(result)
}