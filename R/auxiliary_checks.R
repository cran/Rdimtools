# LIST OF CHECKER FUNCTIONS
# 01. check_ndim   : universal
# 02. check_WMAT   : REE





# 01. check_ndim ----------------------------------------------------------
#' @noRd
#' @keywords internal
check_ndim <- function(ndim, p){
  if ((length(ndim)!=1)||(!is.numeric(ndim))||(ndim<1)||(ndim>p)||is.infinite(ndim)||is.na(ndim)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}


# 02. check_WMAT ----------------------------------------------------------
#' @noRd
#' @keywords internal
check_WMAT <- function(W, n){
  # 1. size argument
  cond1 = ((is.matrix(W))&&(nrow(W)==n)&&(ncol(W)==n))
  # 2. no negative values
  cond2 = (all(W>=0))
  # 3. no Inf of NA
  cond3 = ((!any(is.na(W)))&&(!any(is.infinite(W))))

  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
