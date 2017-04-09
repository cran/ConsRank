#' Branch-and-Bound algorithm to find the median ranking in the space of full (or complete) rankings.
#'
#' Branch-and-bound algorithm to find consensus ranking as definned by D'Ambrosio et al. (2015). If the number of objects to be ranked is large (greater than 20 or 25), it can work for very long time. Use either QuickCons or FASTcons with the option FULL=TRUE instead
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. The data matrix can contain both full and tied rankings, or incomplete rankings. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param PS If PS=TRUE, on the screen some information about how many branches are processed are displayed
#'
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}
#' 
#' @details If the objects to be ranked is large (>25 - 30), it can take long time to find the solutions
#' 
#' @references D'Ambrosio, A., Amodio, S., and Iorio, C. (2015). Two algorithms for finding optimal solutions of the Kemeny rank aggregation problem for full rankings. Electronic Journal of Applied Statistical Analysis, 8(2), 198-213.
#' 
#' @seealso \code{\link{FASTcons}} FAST algorithm algorithm.
#' @seealso \code{\link{QuickCons}} Quick algorithm.
#' 
#' @keywords Median ranking
#' 
#' @examples 
#' data(APAFULL)
#' CR=BBFULL(APAFULL)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export

BBFULL = function(X,Wk=NULL,PS=TRUE)  {
  #Branch and Bound algorithm to find median ranking in the space of full rankings
  #X is a data matrix in which the rows are the judges and the columns indicates the objects
  #Wk is the vector of weigths
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
    
    if (sum(cij==0)==nrow(cij)^2){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    R=findconsensusBB(cij)
    cons1=BBconsensus(R,cij,FULL=TRUE)
    consensus1=cons1$cons
    Po=cons1$pen
    consensus=BBconsensus2(consensus1,cij,Po,PS,FULL=TRUE)
  }
  
  
  if (nrow(consensus)==1) {
    
    Sij=scorematrix(consensus)
    
    if (!is.null(Wk)){
      TauX=sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      TauX=sum(cij*Sij) / (  M*(N*(N-1)) )
    }
    
  } else {
    
    TauX=matrix(0,nrow(consensus),1)
    
    for (k in 1:nrow(consensus)) {
      
      Sij=scorematrix(t(matrix(consensus[k,])))
      
      if (!is.null(Wk)) {
        
        TauX[k,1] = sum(cij*Sij) / ( sum(Wk)*(N*(N-1)) )
        
      } else {
        
        TauX[k,1] = sum(cij*Sij) / (M*(N*(N-1)))
        
      }
      
    }
    
  }
  toc = proc.time()[3]
  colnames(consensus)=colnames(X) 
  #consensus=reordering(consensus)
  eltime=toc-tic
  return(list(Consensus=reordering(consensus), Tau=TauX, Eltime=eltime) )
}
