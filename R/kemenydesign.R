#' Auxiliary function
#'
#' Define a design matrix to compute Kemeny distance
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects represented by the columns.
#' 
#' @return Design matrix
#'
#' @references D'Ambrosio, A. (2008). Tree based methods for data editing and preference rankings. Unpublished PhD Thesis. Universita' degli Studi di Napoli Federico II.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export
#' 
#' @import gtools

kemenydesign <- function(X) {
  
  #Compute the design matrix to compute Kemeny distance
  
  if (is(X,"numeric") & !is(X,"matrix")) {
    X<-matrix(X,ncol=length(X))
  }
  
  
  
  M <- ncol(X)
  N <- nrow(X)
  indice<-combinations(M,2)
  KX<-mat.or.vec(N,(M*(M-1)/2) )
  for (j in 1:nrow(indice)) {
    KX[,j]<-sign(X[,indice[j,1]] - X[,indice[j,2]])*-1
  }
  KX
}
