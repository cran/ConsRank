#'Branch-and-bound algorithm to find consensus (median) ranking according to the Kemeny's axiomatic approach
#'
#'Branch-and-bound algorithm to find consensus ranking as definned by Emond and Mason (2002). If the number of objects to be ranked is large (greater than 15 or 20, specially if there are missing rankings), it can work for very long time.
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param Wk Optional: the frequency of each ranking in the data
#' @param PS If PS=TRUE, on the screen some information about how many branches are processed are displayed
#' 
#' @details This function is deprecated and it will be removed in the 
#' next release of the package. Use function 'consrank' instead.
#' 
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
#' @seealso \code{\link{consrank}} 
#'
#' @keywords Consensus ranking
#' @keywords median ranking
#' 
#' @export




EMCons <- function(X,Wk=NULL,PS=TRUE) {
  .Deprecated(msg = "'EMCons' will be removed in the next release of the package")
  out=consrank(X,wk=Wk,ps=PS,proc=TRUE)
  return(out)
  
}

