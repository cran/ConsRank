#' Kemeny distance
#'
#' Compute the Kemeny distance of a data matrix containing preference rankings, or compute the kemeny distance between two (matrices containing) rankings.
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. If there is only X as input, the output is a square distance matrix
#' @param Y A row vector, or a n by M data matrix in which there are n judges and the same M objects as X to be judged.
#'
#' @return If there is only X as input, d = square distance matrix. If there is also Y as input, d = matrix with N rows and n columns.
#'
#' @references Kemeny, J. G., & Snell, L. J. (1962). Preference ranking: an axiomatic approach. Mathematical models in the social sciences, 9-23.
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @seealso \code{\link{tau_x}} TauX rank correlation coefficient
#' @seealso \code{\link{iw_kemenyd}} item-weighted Kemeny distance
#' 
#' @examples
#' data(Idea)
#' RevIdea<-6-Idea ##as 5 means "most associated", it is necessary compute the reverse 
#' #ranking of each rankings to have rank 1 = "most associated" and rank 5 = "least associated"
#' KD<-kemenyd(RevIdea)
#' KD2<-kemenyd(RevIdea[1:10,],RevIdea[55,])
#' 
#' @keywords Kemeny distance
#' 
#' @export
#' 
#' @import proxy


kemenyd <- function(X,Y=NULL) {
  
  ##Kemeny Distance
  
  if (is(X,"numeric") & !is(X,"matrix")) {
    X<-matrix(X,ncol=length(X))
  }
  
  if (is(Y,"NULL")) {
    X <- kemenydesign(X)
    d<-dist(X,"manhattan")
  } else {
    
    if (is(Y,"numeric") & !is(Y,"matrix")) {
      Y<-matrix(Y,ncol=length(Y))
    }
    
    
    X<-kemenydesign(X)
    Y<-kemenydesign(Y)
    d<-dist(X,Y,"manhattan")
  }
  d
}
