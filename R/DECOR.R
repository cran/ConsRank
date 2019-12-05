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
#' @details This function is deprecated and it will be removed in the 
#' next release of the package. Use function 'consrank' instead.
#' 
#' @references D'Ambrosio, A., Mazzeo, G., Iorio, C., and Siciliano, R. (2017). A differential evolution algorithm for finding the median ranking under the Kemeny axiomatic approach. Computers and Operations Research, vol. 82, pp. 126-138. 
#' 
#' @keywords Differential evolution
#' @keywords Median ranking
#' @keywords Genetic algorithms
#' 
#' @seealso \code{\link{consrank}}
#' 
#' @examples 
#' #not run
#' #data(EMD)
#' #CR=DECOR(EMD[,1:15],EMD[,16])
#' 
#' @export


DECOR <- function(X,Wk=NULL,NP=15,L=100,FF=0.4,CR=0.9,FULL=FALSE) {
  .Deprecated(msg = "'FASTDECOR' will be removed in the next release of the package")
  out=consrank(X,wk=Wk,algorithm="decor",full=FULL,itermax=1,np=NP,gl=L,ff=FF,cr=CR)   
  return(out)
  
}