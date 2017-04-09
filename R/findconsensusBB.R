#' Find a first ranking candidate to be the median ranking
#'
#' Auxiliary function: it finds a first ranking candidate to be the consensus ranking by observing the Combined input matrix
#'
#' @param cij combined input matrix
#'
#' @return candidate median ranking
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export
#' 
#' @import gtools


findconsensusBB = function(cij) {
  
  X=mat.or.vec(1,ncol(cij)) + 1
  N=ncol(X)
  indici=combinations(N,2)
  for (j in 1:nrow(indici)) {
    if ( sign(cij[indici[j,1],indici[j,2]]) == 1 & sign(cij[indici[j,2],indici[j,1]]) == -1 ) {
      X[indici[j,1]]=X[indici[j,1]]+1
    } else if ( sign(cij[indici[j,1],indici[j,2]]) == -1 & sign(cij[indici[j,2],indici[j,1]]) == 1 ) {
      X[indici[j,2]]=X[indici[j,2]] + 1
    } else if (sign(cij[indici[j,1],indici[j,2]]) == -1 & sign(cij[indici[j,2],indici[j,1]]) == -1 ) {
      X[indici[j,1]]= NA
    } else if (sign(cij[indici[j,1],indici[j,2]]) == 1 & sign(cij[indici[j,2],indici[j,1]]) == 1 ){
      X[indici[j,1]]=X[indici[j,1]]+1
      X[indici[j,2]]=X[indici[j,2]] + 1
    }
    
  }
  
  X=(N+1)-X;
  return(X)
}
