#' Auxiliary function
#'
#' Define a design matrix to compute Kemeny distance
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects represented by the columns.
#'
#' @return Design matrix
#'
#' @references D'Ambrosio, A. (2008). Tree based methods for data editing and preference rankings. Unpublished PhD Thesis. Universita' degli Studi di Napoli Federico II.
#'
#' @author Antonio D'Ambrosio \email{antdambr@@unina.it}
#'
#' @details This function now uses an optimized C++ implementation (\code{\link{kemenydesign_cpp}})
#'   for improved performance. The original R implementation is preserved for reference.
#'
#' @seealso \code{\link{kemenydesign_cpp}}
#'
#' @export

kemenydesign <- function(X) {
  # Use C++ implementation for better performance
  return(kemenydesign_cpp(X))
}

# Original R implementation (preserved for reference)
kemenydesign_r <- function(X) {
  if (is(X, "numeric") & !is(X, "matrix")) {
    X <- matrix(X, ncol = length(X))
  }

  KX = sign(X[, 1] - X[, -1]) * -1
  M <- ncol(X)
  for (j in 2:(M - 1)) {
    if (is(KX, "matrix")) {
      KX = cbind(KX, sign(X[, j] - X[, -c(1:j)]) * -1)
    } else {
      KX <- c(KX, sign(X[, j] - X[, -c(1:j)]) * -1)
    }
  }
  if (is(KX, "numeric") & !is(KX, "matrix")) {
    KX <- matrix(KX, ncol = length(KX))
  }

  colnames(KX) <- NULL

  KX
}
