#' Quick algorithm to find up to 4 solutions to the consensus ranking problem
#'
#' The Quick algorithm finds up to 4 solutions. Solutions reached are most of the time optimal solutions. 
#'
#' @param X A N by M data matrix in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once in the sample. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings.
#' @param PS Default PS=FALSE. If PS=TRUE the number of evaluated branches is diplayed
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
#' @examples
#' data(EMD)
#' CR=QuickCons(EMD[,1:15],EMD[,16])
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Amodio, S., D'Ambrosio, A. and Siciliano, R. (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, 249(2), 667-676.
#' 
#' @seealso \code{\link{consrank}}
#'
#' @keywords Quick algorithm
#' 
#' @export

QuickCons <- function(X,Wk=NULL, FULL=FALSE,PS=FALSE) {
  .Deprecated(msg = "'QuickCons' will be removed in the next release of the package")
  out=consrank(X,wk=Wk,algorithm="quick",ps=PS,full=FULL)
  return(out)
  
}