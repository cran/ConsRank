#' Combined input matrix of a data set
#' 
#' Compute the Combined input matrix of a data set as defined by Emond and Mason (2002)
#' 
#' @param X A data matrix N by M, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' 
#' @return The M by M combined input matrix
#' 
#' @examples
#' data(APAred) 
#' CI<-combinpmatr(APAred) 
#' TR<-tabulaterows(APAred) 
#' CI<-combinpmatr(TR$X,TR$Wk)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' 
#' @seealso \code{\link{tabulaterows}} frequency distribution of a ranking data.
#' 
#' @export


combinpmatr <- function (X,Wk=NULL) {
  
  ### COMBINED INPUT MATRIX as defined by Emond and Mason
  
  if (is(Wk,"NULL")) {
    #X must be data matrix with n judges (on the rows) ranking m objects (on the columns)
    CI=matrix(0,ncol(X), ncol(X))
    colnames(CI)=colnames(X)
    row.names(CI)=colnames(X)
    for (i in 1:nrow(X)){
      sm=scorematrix(t(as.matrix(X[i,])))
      CI=CI+sm
    }
  } else {
    
    if (is(Wk,"numeric")) {
      Wk=matrix(Wk,ncol=1)
    }
    
    CI=matrix(0,ncol(X), ncol(X))
    colnames(CI)=colnames(X)
    row.names(CI)=colnames(X)    
    for (i in 1:nrow(X)){
      sm=scorematrix(t(as.matrix(X[i,])))*Wk[i]
      CI=CI+sm
    }
  }
  CI
}
