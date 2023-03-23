#' Branch-and-Bound algorithm to find the median ranking in the space of full (or complete) rankings.
#'
#' Branch-and-bound algorithm to find consensus ranking as defined by D'Ambrosio et al. (2015). If the number of objects to be ranked is large (greater than 20 or 25), it can work for very long time. Use either QuickCons or FASTcons with the option FULL=TRUE instead
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. The data matrix can contain both full and tied rankings, or incomplete rankings. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param PS If PS=TRUE, on the screen some information about how many branches are processed are displayed
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
#' @details If the objects to be ranked is large (>25 - 30), it can take long time to find the solutions
#' 
#' @references D'Ambrosio, A., Amodio, S., and Iorio, C. (2015). Two algorithms for finding optimal solutions of the Kemeny rank aggregation problem for full rankings. Electronic Journal of Applied Statistical Analysis, 8(2), 198-213.
#' 
#' @seealso \code{\link{consrank}}
#' 
#' @keywords Median ranking
#' 
#' @examples 
#' #data(APAFULL)
#' #CR=BBFULL(APAFULL)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export


BBFULL <- function(X,Wk=NULL,PS=TRUE) {
  .Deprecated(msg = "'BBFULL' will be removed in the next release of the package")
  out=consrank(X,wk=Wk,ps=PS,proc=TRUE,full=TRUE)
  return(out)
  
}