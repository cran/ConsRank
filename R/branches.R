#' Auxliary code recalled by other routines
#'
#' Branches discovery 
#'
#' @param brR Current processed branche of the BB algorithm
#' @param cij Combined input matrix
#' @param b Other inputs recalled by main functions
#' @param Po Other inputs recalled by main functions
#' @param ord Other inputs recalled by main functions
#' @param Pb Other inputs recalled by main functions
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. TRUE=TRUE if the function is called by BBFULL algorithm.
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' cR \tab  \tab ranking belonging to the branche\cr
#' pcR \tab       \tab penalty of the current ranking}
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export


branches = function(brR,cij,b,Po,ord,Pb,FULL=FALSE) {
  
  candidate = findbranches(brR,ord,b,FULL)
  Pb =matrix( rep(Pb,nrow(candidate)))
  
  CR=mat.or.vec(nrow(candidate),ncol(candidate))
  addpenalty=matrix(0,nrow(candidate),1)
  QR=mat.or.vec(nrow(candidate),ncol(candidate))
  
  for (gm in 1:nrow(candidate)) {
    
    CR[gm,]=candidate[gm,]
    addpenalty[gm,]=PenaltyBB2(cij,candidate[gm,],ord[b])
    
    if (Pb[gm]+addpenalty[gm] > Po) {
      
      CR[gm,]=-10.0e+15
      addpenalty[gm]=-10.0e+15
      
    }
    QR[gm,]=CR[gm,]
  }
  Pbr=addpenalty+Pb
  idp=Pbr<0
  
  if (sum(idp)==0) {
    
    R=QR
    
  } else if (sum(idp==F)==nrow(QR)) {
    
    Pbr=NULL
    Pb=NULL
    R=NULL
    
  } else {
    Pbr=t(matrix(Pbr[idp==FALSE,],1))
    if (sum(idp==F)==1) {
      R=t(matrix(QR[idp==FALSE,]))
    } else {
      R=QR[idp==FALSE,]
    }
  }
  
  return(list(cR=R,pcR=Pbr))
}
