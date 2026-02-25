#' Score matrix according Emond and Mason (2002)
#'
#' Given a ranking, it computes the score matrix as defined by Emond and Mason (2002)
#'
#' @param X a ranking (must be a row vector or, better, a matrix with one row and M columns)
#' @param use_cpp Logical. If TRUE (default), uses C++ implementation for faster computation. If FALSE, uses the R implementation.
#'
#' @return the M by M score matrix
#'
#' @examples
#' Y <- matrix(c(1,3,5,4,2),1,5)
#' SM<-scorematrix(Y)
#' #
#' Z<-c(1,2,4,3)
#' SM2<-scorematrix(Z)
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' 
#' @seealso \code{\link{combinpmatr}} The combined inut matrix
#' 
#' @export


scorematrix <- function (X, use_cpp=TRUE) {
  
  # Validate / standardize input
  itemnames <- if (is.matrix(X)) colnames(X) else names(X)
  
  if (is.matrix(X)) {
    if (nrow(X) != 1L) stop("X must be a row vector (1 x m) or a numeric vector")
    x <- as.numeric(X[1, , drop = TRUE])
  } else if (is.numeric(X)) {
    x <- as.numeric(X)
  } else {
    stop("X must be a numeric vector or a 1-row numeric matrix")
  }
  
  if (use_cpp) {
    # Chiama scorematrix_cpp() - C++ ottimizzato
    sm <- scorematrix_cpp(matrix(x, nrow = 1))
    # Restore names...
    return(sm)
  }
  
  # Fallback R implementation
  

  
  m <- length(x)
  if (m == 0L) return(matrix(numeric(0), 0, 0))
  
  # Vectorized: build logical matrix of comparisons x_i <= x_j
  # NA entries in 'leq' correspond to comparisons with NA in x
  leq <- outer(x, x, "<=")
  
  # Initialize result matrix with -1 and set 1 where x_i <= x_j
  sm <- matrix(-1L, nrow = m, ncol = m)
  sm[leq] <- 1L
  
  # Diagonal to 0
  diag(sm) <- 0L
  
  # Where comparison is NA (e.g., x contains NA), set to 0
  sm[is.na(leq)] <- 0L
  
  # Restore names if present
  if (!is.null(itemnames) && length(itemnames) == m) {
    rownames(sm) <- itemnames
    colnames(sm) <- itemnames
  }
  
  sm
}