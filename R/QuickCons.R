#' Quick algorithm to find up to 4 solutions to the consensus ranking problem
#'
#' The Quick algorithm finds up to 4 solutions. Solutions reached are most of the time optimal solutions. 
#'
#' @param X A N by M data matrix in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once in the sample. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings.
#' @param PS Default PS=FALSE. If PS=TRUE the number of evaluated branches is diplayed
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}
#' 
#' @examples
#' data(EMD)
#' CR=QuickCons(EMD[,1:15],EMD[,16])
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Sonia Amodio \email{sonia.amodio@unina.it}
#' 
#' @references Amodio, S., D'Ambrosio, A. and Siciliano, R. (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, 249(2), 667-676.
#' 
#' @seealso \code{\link{FASTcons}} FAST algorithm.
#' @seealso \code{\link{QuickCons}} Quick algorithm.
#'
#' @keywords Quick algorithm
#' 
#' @export




QuickCons = function(X,Wk=NULL, FULL=FALSE,PS=FALSE)   {
  
  
  if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
  }
  
  M = nrow(X)
  N=ncol(X)
  tic = proc.time()[3]
  
  if (M==1) {
    consensus = X
    TauX = 1
  } else {
    if (!is.null(Wk)) {
      
      if (is.numeric(Wk)) {
        Wk=matrix(Wk,ncol=1)
      }
      
      cij = combinpmatr(X,Wk)
    } else {
      cij = combinpmatr(X)
    }
    
    if (sum(cij==0)==length(cij)){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    if (sum(sign(cij+diag(N)))==length(cij)){
      print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
      return()
      
    }
    
    R=findconsensusBB(cij)
    R1=(N+1)-R
    consensusA = BBconsensus(R,cij,FULL,PS)$cons
    consensusB = BBconsensus(consensusA,cij,FULL,PS)$cons
    consensusC = BBconsensus(R1,cij,FULL,PS)$cons
    consensusD = BBconsensus(consensusC,cij,FULL,PS)$cons
    consensus = unique(reordering(rbind(consensusA,consensusB,consensusC,consensusD)))
    howcons = nrow(consensus)
    
    
  }
  #d=kemenyd(X,consensus$cons)
  
  Taux=matrix(0,nrow(consensus),1)
  for (k in 1:nrow(consensus)) {
    
    #Sij=scorematrix(t(as.matrix(consensus[k,])))
    Sij=scorematrix(matrix(consensus[k,],1,N))
    
    if (!is.null(Wk)){
      Taux[k,1]=sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      Taux[k,1]=sum(cij*Sij) / (  M*(N*(N-1)) )
    }
    
  }
  
  if (howcons>1) {
    nco=which(Taux==max(Taux))
    if (length(nco)>1) {
      consensus=consensus[nco,]
      Taux=matrix(rep(max(Taux),nrow(consensus)),nrow(consensus),1)
    } else {
      Taux=max(Taux)
      #consensus = t(matrix(consensus[nco,]))
      consensus = matrix(consensus[nco,],1,N)
    }
  }
  colnames(consensus)=colnames(X) 
  toc = proc.time()[3]
  eltime=toc-tic
  return(list(Consensus=reordering(consensus), Tau=Taux, Eltime=eltime) )
}
