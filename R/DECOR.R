#' Differential Evolution algorithm for Median Ranking
#'
#' Differential evolution algorithm for median ranking detection. 
#' It works with full, tied and partial rankings. The solution con be
#' constrained to be a full ranking or a tied ranking
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param NP The number of population individuals
#' @param L Generations limit: maximum number of consecutive generations without improvement
#' @param FF The scaling rate for mutation. Must be in [0,1]
#' @param CR The crossover range. Must be in [0,1]
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings.
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo \email{giuliomazzeo@gmail.com}
#' 
#' @details It works with a very large number of items to be ranked. Empirically, the number
#' of population individuals (the NP parameter) can be set equal to 10, 20 or 30
#' for problems till 20, 50 and 100 items. Both scaling rate and vrossover rati (parameters
#' FF and CR) must be set by the user. The default options (FF=0.4, CR=0.9) work well
#' for a large variety of data sets
#' 
#' @references D'Ambrosio, A., Mazzeo, G., Iorio, C., and Siciliano, R. (2017). A differential evolution algorithm for finding the median ranking under the Kemeny axiomatic approach. Computers and Operations Research, vol. 82, pp. 126-138. 
#' 
#' @keywords Differential evolution
#' @keywords Median ranking
#' @keywords Genetic algorithms
#' 
#' @seealso \code{\link{FASTcons}} FAST algorithm.
#' @seealso \code{\link{QuickCons}} Quick algorithm.
#' @seealso \code{\link{EMCons}} Branch-and-bound algorithm.
#' 
#' @examples 
#' data(EMD)
#' CR=DECOR(EMD[,1:15],EMD[,16])
#' 
#' @export
#' 


DECOR = function(X,Wk=NULL,NP=15,L=100,FF=0.4,CR=0.9,FULL=FALSE){

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
    
    
    COR=DECORcore(cij,NJ,NP,L,FF,CR,FULL)
    
    Consensus=COR$ConsR
    TauX=COR$Tau
    
    
    
    
  }
  toc=proc.time()[3]
  eltime=toc-tic
  return(list(Consensus=reordering(Consensus), Tau=TauX, Eltime=eltime) )
}
