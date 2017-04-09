#' Auxiliary function
#'
#' Assign a penalty to the branches of the BB algorithms
#'
#' @param cij combined input matrix
#' @param candidate candidate to be the median ranking
#' @param ord other input called by other functions
#'
#' @return computed penalty
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export

PenaltyBB2 = function(cij,candidate,ord)   #indice must be order(CR)
  
  ## DETERMINATION OF PENALTIES FOR THE BRANCH AND BOUND ALGORITHM
  
{
  
  Ds=t(mat.or.vec(1,(length(ord)-1)));
  addpenalty=t(mat.or.vec(1,(length(ord)-1)));
  
  for (k in 1:(length(ord)-1)) {
    
    Ds[k,1]=sign(candidate[ord[length(ord)]]-candidate[ord[k]]);
    
    if (Ds[k,1]==1) {
      
      
      if ( sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
        addpenalty[k,1]=cij[ord[length(ord)],ord[k]]-cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 & sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 & sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]=cij[ord[length(ord)],ord[k]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 & sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]=0
      }
    }
    else if (Ds[k,1]==-1) {
      if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
        addpenalty[k,1]=0
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]=cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]=cij[ord[k],ord[length(ord)]]-cij[ord[length(ord)],ord[k]]
      }
    }
    
    else if (Ds[k,1]==0) {
      if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
        addpenalty[k,1]=-cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]=0
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]=-cij[ord[length(ord)],ord[k]]
      }
    }
  }
  
  addpenalty=sum(addpenalty)
  
}
