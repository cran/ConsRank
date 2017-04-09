#' Auxiliary function
#'
#' Find correct branches in the Branch-and-Bound algorithms
#'
#' @param R Candidate to be the consensus ranking
#' @param ord other input values recalled by other routines
#' @param b other input values recalled by other routines
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. TRUE=TRUE if the function is called by BBFULL algorithm. 
#'
#' @return a candidate to be the median ranking
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Sonia Amodio \email{sonia.amodio@unina.it}
#' 
#' @export


findbranches = function(R,ord,b,FULL=FALSE) {
  
  KR=t(matrix(R[ord[b]]))
  KR=KR[-length(KR)]
  MO=max(KR)
  MI=min(KR)
  aa=1
  KO=1
  KR[length(KR)+1]=MO+1;
  R[ord[b]]=KR
  candidate=mat.or.vec(nrow(R), ncol(R))
  
  while (KO==1)  {
    candidate=rbind(candidate,R)
    
    if (aa==1){
      candidate=matrix(candidate[-1,],1,ncol(candidate))
    }
    if (FULL==FALSE){
      R[ord[b[length(b)]]]=R[ord[b[length(b)]]]-1 }else{
        R[ord[b[length(b)]]]=R[ord[b[length(b)]]]-2
      }
    
    if (MI-R[ord[b[length(b)]]] > 1) {
      
      KO=0
      
    }
    
    aa=aa+1
    
  }
  
  Rt=candidate
  
}
