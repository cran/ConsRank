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
#' @details This function is deprecated and it will be removed in the 
#' next release of the package. Use function 'consrank' instead.
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
#' @seealso \code{\link{consrank}}
#' 
#' @examples 
#' #data(EMD)
#' #CR=FASTDECOR(EMD[,1:15],EMD[,16])
#' 
#' @keywords Median ranking
#' @keywords Differential evolution
#' 
#' @export


FASTDECOR <- function(X,Wk=NULL,maxiter=10,NP=15,L=100,FF=0.4,CR=0.9,FULL=FALSE,PS=TRUE) {
  .Deprecated(msg = "'FASTDECOR' will be removed in the next release of the package")
  out=consrank(X,wk=Wk,ps=PS,algorithm="decor",full=FULL,itermax=maxiter,np=NP,gl=L,ff=FF,cr=CR)   
  return(out)
  
}