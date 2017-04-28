#'Branch-and-bound algorithm to find consensus (median) ranking according to the Kemeny's axiomatic approach
#'
#'Branch-and-bound algorithm to find consensus ranking as definned by Emond and Mason (2002). If the number of objects to be ranked is large (greater than 15 or 20, specially if there are missing rankings), it can work for very long time.
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param PS If PS=TRUE, on the screen some information about how many branches are processed are displayed
#' 
#' @details If the objects to be ranked is large (>15-20) with some missing, it can take long time to find the solutions. If the searching space is limited to the space of full rankings (also incomplete rankings, but without ties), use the function BBFULL or the functions FASTcons and QuickCons with the option FULL=TRUE.
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}
#' 
#' @examples
#' data(Idea)
#' RevIdea=6-Idea 
#' # as 5 means "most associated", it is necessary compute the reverse ranking of 
#' # each rankings to have rank 1 = "most associated" and rank 5 = "least associated"
#' CR=EMCons(RevIdea)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' 
#' @seealso \code{\link{FASTcons}} FAST algorithm algorithm.
#' @seealso \code{\link{QuickCons}} Quick algorithm.
#' @seealso \code{\link{BBFULL}} Branc-and-bound algorithm for full rankings.
#'
#' @keywords Consensus ranking
#' @keywords median ranking
#' 
#' @export



EMCons = function(X,Wk=NULL,PS=TRUE)  {
  #Emond and Mason Branch and Bound algorithm to find median ranking
  #X is a data matrix in which the rows are the judges and the columns indicates the objects
  #Wk is the vector of weigths
  if (class(X)=="data.frame") {
    #colnames(X)=NULL
    X=as.matrix(X)
  }
  
  
  
  M = nrow(X)
  N=ncol(X)
  callps=PS
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
    cons1=BBconsensus(R,cij,FULL=FALSE,PS=FALSE)
    consensus1=cons1$cons
    Po=cons1$pen
    consensus=BBconsensus2(consensus1,cij,Po,PS=callps,FULL=FALSE)
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

