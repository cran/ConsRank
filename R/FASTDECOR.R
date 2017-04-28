#' FAST algorithm calling DECOR
#'
#' FAST algorithm repeats DECOR a prespecified number of time. 
#' It returns the best solutions among the iterations
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param maxiter maximum number of iterations. Default 10
#' @param NP The number of population individuals
#' @param L Generations limit: maximum number of consecutive generations without improvement
#' @param FF The scaling rate for mutation. Must be in [0,1]
#' @param CR The crossover range. Must be in [0,1]
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings. In this case, the data matrix must contain full rankings.
#' @param PS Default PS=TRUE. If PS=TRUE the number of a multiple of 5 iterations is diplayed
#'
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}
#' 
#' @references D'Ambrosio, A., Mazzeo, G., Iorio, C., and Siciliano, R. (2017). A differential evolution algorithm for finding the median ranking under the Kemeny axiomatic approach. Computers and Operations Research, vol. 82, pp. 126-138. 
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo \email{giuliomazzeo@gmail.com}
#' 
#' @examples 
#' #data(EMD)
#' #CR=FASTDECOR(EMD[,1:15],EMD[,16])
#' 
#' @keywords Median ranking
#' @keywords Differential evolution
#' 
#' @export


FASTDECOR = function(X,Wk=NULL,maxiter=10,NP=15,L=100,FF=0.4,CR=0.9,FULL=FALSE,PS=TRUE){
  
  #check if X is a matrix
  if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
  }
  
  M = nrow(X)
  N=ncol(X)
  tic = proc.time()[3]  
  
  #check if there are trivial solutions
  
  if (M==1) { 
    consensus = X
    TauX = 1
    
  } else {
    
    if (!is.null(Wk)) {
      
      if (is.numeric(Wk)) {
        Wk=matrix(Wk,ncol=1)
        NJ=sum(Wk)
      }
      
      cij = combinpmatr(X,Wk)
      
    } else {
      
      cij = combinpmatr(X)
      NJ=nrow(X)
    }
    
    
    if (sum(cij==0)==length(cij)){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    if (sum(sign(cij+diag(N)))==length(cij)){
      print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
      return()
      
    }
    
    
    sol=matrix(0,1,N)
    taos=0
    
    for (iter in 1:maxiter){
    COR=DECORcore(cij,NJ,NP,L,FF,CR,FULL)
    
    sol=rbind(sol,COR$ConsR)
    taos=rbind(taos,COR$Tau)
    
    if (PS==TRUE) {
      
      if (iter%%5==0){
      
      dsp1=paste("Iteration",iter,sep=" ")
      print(dsp1)
      }
      
    }

    
    }
    
  }
  
  sol=sol[-1,]
  taos=matrix(taos[-1],(length(taos)-1),1)
  
  if (is.null(nrow(sol))){sol=matrix(sol,1,N)}
  
  bestindex=which(taos==max(taos))
  sol=sol[bestindex,]
  taos=max(taos)
  
  if (is.null(nrow(sol))){
    sol=matrix(sol,1,N)
    Consensus=sol
  } else {
    Consensus=unique(sol)
  }
  
  if (is.null(nrow(Consensus))){
    Consensus=matrix(Consensus,1,N)
    }

  TauX=matrix(rep(taos),nrow(Consensus),1)
  colnames(Consensus)=colnames(X) 
  row.names(Consensus)=NULL
  
  
  toc=proc.time()[3]
  eltime=toc-tic
  return(list(Consensus=reordering(Consensus), Tau=TauX, Eltime=eltime) )
}
