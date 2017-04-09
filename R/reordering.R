#' Given a vector (or a matrix), returns an ordered vector (or a matrix with ordered vectors)
#'
#' Given a ranking of M objects (or a matrix with M columns), it reduces it in "natural" form (i.e., with integers from 1 to M)
#'
#' @param X a ranking, or a ranking data matrix
#'
#' @return a ranking in natural form
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export



reordering = function (X) {
  
  if (nrow(X)==1) {
    G=X
    OX = order(X)
    SX=matrix(sort(X),nrow=1)
    SX=SX-min(SX)+1
    DC=rbind(0,diff(t(SX)))
    
    
    for (i in 1:(ncol(X)-1)) {
      if (DC[i+1,1] >= 1) {
        SX[1,i+1]=SX[1,i]+1
      } else if (DC[i+1,1] == 0) {
        SX[1,i+1]=SX[1,i]
      }
    }
    
    G[1,OX]=SX
    
  } else {
    
    G=X
    for (j in 1:nrow(X)) {
      
      OX = order(X[j,])
      SX=matrix(sort(X[j,]),nrow=1)
      SX=SX-min(SX)+1
      DC=rbind(0,diff(t(SX)))
      
      for (i in 1:(ncol(X)-1)) {
        
        if (DC[i+1,1] >= 1) {
          SX[1,i+1]=SX[1,i]+1
        } else if (DC[i+1,1] == 0) {
          SX[1,i+1]=SX[1,i]
        }
      }
      G[j,OX]=SX
    }
    
  }
  
  G
}
