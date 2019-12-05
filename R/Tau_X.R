#' TauX (tau exstension) rank correlation coefficient
#'
#' Tau exstension is a new rank correlation coefficient defined by Emond and Mason (2002)
#'
#' @param X a M by N data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. If there is only X as input, the output is a square matrix containing the Tau_X rcc.
#' @param Y A row vector, or a n by M data matrix in which there are n judges and the same M objects as X to be judged.
#'
#' @return Tau_x rank correlation coefficient
#'
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @seealso \code{\link{kemenyd}} Kemeny distance
#' 
#' @keywords TauX rank correlation coefficient
#' 
#' @examples
#' data(BU)
#' RD<-BU[,1:3]
#' Tau<-tau_x(RD)
#' Tau1_3<-tau_x(RD[1,],RD[3,])
#' 
#' @export 


tau_x <- function(X,Y=NULL) {
  
  if(is(X,"numeric") & !is(X,"matrix")){
    X<-matrix(X,ncol=length(X))
  }
  
  
  N<-ncol(X)
  maxd<-N*(N-1)
  if (is(Y,"NULL")){
    d<-kemenyd(X)
    tau <- 1-(2*d/maxd)
  } else {
    if(is(Y,"numeric") & !is(Y,"matrix")){
      Y<-matrix(Y,ncol=length(Y))
    }
    d<-kemenyd(X,Y)
    tau <- 1-(2*d/maxd)
  }
  tau 
}

#' @rdname tau_x
#' @export
Tau_X <- tau_x

