#' Find the first approximation to the consensus ranking. Most of the time the output is a solution, maybe not unique
#' 
#' Find a first approximation to the consensus ranking.
#' 
#' @param RR Candidate to be the consensus ranking
#' @param cij Combined input matrix of the data set
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. FULL=TRUE if the function is called by BBFULL algorithm.
#' @param PS Default PS=FALSE. If PS=TRUE the number of evaluated branches is diplayed
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' cons \tab  \tab a first approximation of the median ranking\cr
#' pen \tab       \tab penalty value}
#' 
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Amodio, S., D'Ambrosio, A. and Siciliano, R. (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, 249(2), 667-676.
#' 
#' @export


BBconsensus = function(RR,cij,FULL=FALSE,PS=FALSE) {
  
  #Branch and Bound Algorithm to find the the consensus ranking PART I As modified by D'AMBROSIO (2008).
  #Find the first approximation to the consensus ranking. Most of the time CR
  #is a solution, maybe not unique
  #Input:
  #       RR -> First solution candidate to be the consensus ranking
  #       cij -> Combined Input Matrix of the M individuals which judge n
  #       objects
  #Output:
  #       CR -> Consensus Ranking
  #
  #
  #References: Amodio et al.,2015; D'Ambrosio et al., 2016.
  
  CR=RR
  sij=scorematrix(RR)
  Po = sum(abs(cij))-sum(cij*sij)
  a = t(matrix(sort(RR,decreasing = TRUE)))
  ord = t(matrix(order(RR,decreasing = TRUE)))
  R=RR
  addpenalty=matrix(0,length(a),1)
  
  # exploration of the initial solution
  for (k in 2:length(a)) {
    #print(k)
    b = 1:k
    R = ReorderingBB(R)
    KR=t(matrix(R[ord[b]]))
    KR=KR[-length(KR)]
    MO=max(KR)
    MI=min(KR)
    aa=1
    KO=1
    KR[length(KR)+1]=MO+1
    R[ord[b]]=KR
    candidate=matrix(0,nrow(RR), ncol(RR))
    Pb = matrix(0, 1, 1)
    while (KO==1)  {
      #browser()
      candidate=rbind(candidate,R)
      #if (ncol(candidate>ncol(RR))) {
      
      #print(dim(candidate))
      #print(dim(R))
      #print(class(candidate))
      #print(class(R))
      
      # }
      if (aa==1){
        candidate=matrix(candidate[-1,],1,ncol(candidate))
      }
      
      Sij=scorematrix(matrix(candidate[aa,],1,ncol(R)))
      # print(Sij)
      #print(candidate)
      #flush.console()
      Pb=rbind(Pb,sum(abs(cij))-sum(cij*Sij))
      if (aa==1) {
        Pb=matrix(Pb[-1,],1,1)
      }
      # print(Pb)
      if (Pb[aa]==0) {
        
        CR = R
        Po = 0
        Pc = 0
        break
      }
      Pc=1
      if(FULL==TRUE){
        R[ord[b[length(b)]]] = R[ord[b[length(b)]]]-2 }else{
          R[ord[b[length(b)]]] = R[ord[b[length(b)]]]-1
        }
      if (MI-R[ord[b[length(b)]]] > 1) {
        KO = 0
      }
      aa=aa+1
      
    }
    
    if (PS==TRUE) {
      
      dsp2=paste("evaluated",nrow(candidate),"branches",sep=" ")
      print(dsp2)
      
    }
    
    if (Pc == 0) {
      break
    }
    
    minp=min(Pb)
    posp=which(Pb==min(Pb))
    
    if (minp<=Po) {
      Po=minp
      CR=t(matrix(candidate[posp[1],]))
      R=CR
      addpenalty[k,1]=PenaltyBB2(cij,R,ord[b])
    } else {
      R = CR
      addpenalty[k,1]=PenaltyBB2(cij,R,ord[b])
    }
    
    candidate = mat.or.vec(nrow(R), ncol(R))
    Pb = mat.or.vec(1, 1)
    
  }
  
  if (Pc==0) {
    Po=0
    addpenalty = 0
  }  else {
    Poo=sum(addpenalty)
  }
  
  SIJ = scorematrix(CR)
  Po=sum(addpenalty)
  
  return(list(cons=CR,pen=Po))
}
