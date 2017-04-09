#' Auxiliary function
#'
#' It allows to reord the objects to be processed in the BB algorithms
#'
#' @param RR A ranking
#' 
#' @return A reordered ranking
#'
#' @author Sonia Amodio \email{sonia.amodio@unina.it}
#' 
#' @export

ReorderingBB = function(RR) {
  
  RR = RR+1
  R = RR;
  k = ncol(R)
  neword = order(R)
  indexing = mat.or.vec(1, ncol(R)-1)
  for (j in (k-1):1) {
    indexing[j] = R[neword[j+1]]-R[neword[j]]
  }
  
  if (sum(indexing==0)>0) {
    J = 1
    while (J<=ncol(indexing)) {
      if (indexing[J]==0) {
        R[neword[J+1]]=R[neword[J]]
        J=J+1
      }
      else if (indexing[J]>0) {
        R[neword[J+1]] = R[neword[J]]+2
        J=J+1
      }
    }
  }
  else  {
    J = 1
    while (J<= ncol(indexing)) {
      R[neword[J+1]] = R[neword[J]] + 2
      J=J+1
    }
  }
  R
}

#-------------------------------------------------------------------------------
