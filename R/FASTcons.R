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
#' ##data(EMD)
#' ##X=EMD[,1:15]
#' ##Wk=matrix(EMD[,16],nrow=nrow(X))
#' ##CR=FASTcons(X,Wk,maxiter=100)
#' ##These lines produce all the three solutions in less than a minute.
#'
#' data(sports)
#' CR=FASTcons(sports,maxiter=5)
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


FASTcons <- function(X, Wk=NULL, maxiter=50, FULL=FALSE, PS=FALSE) {
  .Deprecated(msg = "'FASTcons' will be removed in the release of the package")
  out=consrank(X,wk=Wk,ps=PS,algorithm="fast",itermax=maxiter)
  return(out)
  
}