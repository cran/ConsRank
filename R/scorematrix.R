#' Score matrix according Emond and Mason (2002)
#'
#' Given a ranking, it computes the score matrix as defined by Emond and Mason (2002)
#'
#' @param X a ranking (must be a row vector or, better, a matrix with one row and M columns)
#'
#' @return the M by M score matrix
#'
#' @examples
#' Y = matrix(c(1,3,5,4,2),1,5)
#' SM=scorematrix(Y)
#' #
#' Z=c(1,2,4,3)
#' SM2=scorematrix(Z)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' 
#' @seealso \code{\link{combinpmatr}} The combined inut matrix
#' 
#' @export

scorematrix = function (X) {
  
  ### SCORE MATRIX OF RANK DATA ACCORDING EMOND AND MASON
  
  if (is.numeric(X) & !is.matrix(X)){
    X=matrix(X,ncol=length(X))
  }
  
  c=ncol(X)
  
  #X must be a row vector containing a ranking of m objects
  sm=matrix(0,c,c)
  
  for (j in 1:c){
    diffs=sign(X[j]-X[setdiff(1:c,j)])
    ind=setdiff(1:c,j)
    sm[j,ind]=diffs
  }
  
  idn=is.na(sm)
  sm=((sm<=0)*2-1)-diag(c)
  sm[idn]=0
  sm
}

