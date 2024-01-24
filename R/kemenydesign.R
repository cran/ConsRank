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


#@import gtools

kemenydesign <- function(X) {
  #to be modified (23/06/2023)
  
  #Compute the design matrix to compute Kemeny distance
  
  if (is(X,"numeric") & !is(X,"matrix")) {
    X<-matrix(X,ncol=length(X))
  }
  
  #--ok
  #nc <- M*(M-1)/2
  #modifica fatta il 16/01/2024, funziona meglio, è più veloce e consente di processare
  #moltissimi oggetti
  KX=sign(X[,1]-X[,-1])*-1
  M <- ncol(X)
  for (j in 2:(M-1)){
    if(is(KX,"matrix")){
      KX=cbind(KX,sign(X[,j]-X[,-c(1:j)])*-1)
    } 
    else 
    {
      KX <- c(KX,sign(X[,j]-X[,-c(1:j)])*-1)
    }
  }
  if(is(KX,"numeric") & !is(KX,"matrix")) KX <- matrix(KX,ncol=length(KX))
  
  colnames(KX) <- NULL
  
  KX
  
  
  #--older code
  # M <- ncol(X)
  # N <- nrow(X)
  # indice<-combinations(M,2)
  # KX<-mat.or.vec(N,(M*(M-1)/2) )
  # for (j in 1:nrow(indice)) {
  #   KX[,j]<-sign(X[,indice[j,1]] - X[,indice[j,2]])*-1
  # }
  # KX
}







