#' Transform a ranking into a ordering. 
#'
#' Given a ranking (or a matrix of rank data), transforms it into an ordering (or a ordering matrix)
#'
#' @param x a  ranking, or a n by m data matrix in which there are n judges ranking m objects
#' @param m the number of objects
#' @param label optional: the name of the objects
#' @param labs labs = 1 displays the names of the objects if there is argument "label", otherwise displays the permutation of first m integer.
#'  labs = 2 is to be used only if the argument "label" is not defined. In such a case it displays the permutation of the first m letters
#'  
#' @details This function is deprecated and it will be removed in the 
#' next release of the package. Use function 'rank2order' instead.
#'
#' @return the ordering
#'
#' @author Sonia Amodio \email{sonia.amodio@unina.it}
#' 
#' @seealso \code{\link{rank2order}}
#' 
#' @examples
#' data(Idea)
#' TR=tabulaterows(Idea)
#' Ord=labels(TR$X,ncol(Idea),colnames(Idea),labs=1)
#' Ord2=labels(TR$X,ncol(Idea),labs=2)
#' cbind(Ord,TR$Wk)
#' cbind(Ord2,TR$Wk)
#' 
#' @export









labels <- function(x, m, label = 1:m, labs ){
  
  ## Place labels in a data matrix X of rankings (N judges by M objects)
  #m is the number of objects
  #label (optional) is the vector of the objects names
  #labs = 1 or 2
  #source('reordering.r')
  # if the class of the object is different from 'matrix' transform it in 'matrix'
  if(!is(x,"matrix")){
    obs <- length(x)
    XX <- matrix(x, ncol = obs)
  } else {
    XX <- x
  }
  
  nj <- nrow(XX)
  nob <- ncol(XX)
  
  ## if length of the object is higher than m, last number is the penalty
  #if(nob > m){
  ## if the number of rows is 1 is a vector
  #  if(nj == 1){
  #    pens = x[m+1]
  #    X = matrix(reordering(XX[1:m]), m, ncol = m)
  #    } else {
  #    pens = x[,m+1]
  #    X = t(apply(x, 1, function(g) reordering(g, m)))
  #    }
  #} else {
  X <- XX
  #}
  if(labs ==1){
    let <- label
  } else if(labs == 2){
    let <- LETTERS[label]
  }
  
  out <- rep(0, nj)
  for(i in 1:nj){
    
    ord <- rank(X[i,])
    orders <- tapply(let, ord, sort)
    
    names1 <- NULL
    for(j in 1:length(orders)){
      if(length(orders[[j]]) > 1){
        nams <- paste('(', paste(orders[[j]], sep = '', collapse = ' '), ')', sep = '', collapse='')
      } else {
        nams <- paste(orders[[j]], collapse = ' ')
      }
      names1 <- c(names1, nams)
    }
    out[i] <- paste(names1, collapse = ' ' )
  }
  out <- matrix(out, nrow = nj)
  
  #if(nob > m){
  #dat = data.frame(data = out, pens = pens)
  #} else {
  dat <- out
  #}
  
  return(dat)
}

