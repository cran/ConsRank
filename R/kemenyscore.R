#' Score matrix according Kemeny (1962)
#'
#' Given a ranking, it computes the score matrix as defined by Emond and Mason (2002)
#'
#' @param X a ranking (must be a row vector or, better, a matrix with one row and M columns)
#'
#' @return the M by M score matrix
#'
#' @examples
#' Y <- matrix(c(1,3,5,4,2),1,5)
#' SM<-kemenyscore(Y)
#' #
#' Z<-c(1,2,3,2)
#' SM2<-kemenyscore(Z)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Kemeny, J and Snell, L. (1962). Mathematical models in the social sciences.
#' 
#' @seealso \code{\link{scorematrix}} The score matrix as defined by Emond and Mason (2002)
#' 
#' @export

kemenyscore <- function (X) {
  
  ### SCORE MATRIX OF RANK DATA ACCORDING TO KEMENY
  
  itemnames<-names(X)
  
  if (is(X,"numeric") & !is(X,"matrix")){
    X<-matrix(X,ncol=length(X))
  }
  
  c<-ncol(X)
  
  #X must be a row vector containing a ranking of m objects
  sm<-matrix(0,c,c)
  colnames(sm)<-itemnames
  row.names(sm)<-itemnames  
  
  for (j in 1:c){
    diffs<-sign(X[j]-X[setdiff(1:c,j)])*-1
    ind<-setdiff(1:c,j)
    sm[j,ind]<-diffs
  }
  
  #sm=((sm<=0)*2-1)-diag(c)
  sm
}

