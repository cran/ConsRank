#' FAST algorithm to find consensus (median) ranking.
#'  
#' FAST algorithm to find consensus (median) ranking defined by Amodio, D'Ambrosio and Siciliano (2016). It returns at least one of the solutions. If there are multiple solutions, sometimes it returns all the solutions, sometimes it returns some solutions, always it returns at least one solution.
#' 
#' @param X is a ranking data matrix
#' @param Wk is a vector of weights
#' @param maxiter maximum number of iterations: default = 50.
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings. 
#' @param PS Default PS=FALSE. If PS=TRUE the number of current iteration is diplayed
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}
#'
#' @examples
#' ##data(EMD)
#' ##X=EMD[,1:15]
#' ##Wk=matrix(EMD[,16],nrow=nrow(X))
#' ##CR=FASTcons(X,Wk,maxiter=100)
#' ##These lines produce all the three solutions in less than a minute.
#'
#' data(sports)
#' CR=FASTcons(sports,maxiter=10)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Sonia Amodio \email{sonia.amodio@unina.it}
#' 
#' @references Amodio, S., D'Ambrosio, A. and Siciliano, R. (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, 249(2), 667-676.
#' 
#' @seealso \code{\link{EMCons}} Emond and Mason branch-and-bound algorithm.
#' @seealso \code{\link{QuickCons}} Quick algorithm.
#'
#' @keywords FAST algorithm
#' 
#' @export


FASTcons = function(X, Wk=NULL, maxiter=50, FULL=FALSE, PS=FALSE)   {
  
  
  if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
  }
  
  M = nrow(X)
  N=ncol(X)
  
  tic = proc.time()[3]
  if (M==1) {
    CR = X
    Taux = 1
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
    
    CR=matrix(0,maxiter,ncol(X))
    for (iter in 1:maxiter) {
      
      if (iter%%2==0) {
        R=matrix(sample(1:ncol(X),replace=T),1,ncol(X))
      } else {
        R=matrix(sample(1:ncol(X)),1,ncol(X))
      }
      
      consensus1 = BBconsensus(R,cij, FULL)
      cons=matrix(consensus1$cons,1,ncol(X))
      consensus = BBconsensus(cons,cij, FULL)
      #print(cons)
      #print(R)
      #flush.console()
      CR[iter,]=matrix(consensus$cons,1,ncol(X))
      if (PS==TRUE) {
        
        dsp1=paste("Iteration",iter,sep=" ")
        print(dsp1)
      }
      
    }
  }
  #d=kemenyd(X,consensus$cons)
  
  Taux=matrix(0,nrow(CR),1)
  for (k in 1:nrow(CR)) {
    Sij=scorematrix(matrix(CR[k,],1,ncol(X)))
    if (!is.null(Wk)){
      Taux[k,]=sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      Taux[k,]=sum(cij*Sij) / (  M*(N*(N-1)) )
    }
  }
  
  CR=reordering(CR)
  indice=which(Taux==max(Taux));
  Taux=max(Taux)
  CR=matrix(CR[indice,],ncol=N)
  if (nrow(CR>1)){
    CR=unique(CR)
  }
  if (!is.null(dim(CR))) {
    Taux=matrix(rep(Taux,nrow(CR)))
  }
  
  colnames(CR)=colnames(X) 
  toc = proc.time()[3]
  eltime=toc-tic
  return(list(Consensus=CR, Tau=Taux, Eltime=eltime) )
}