#' Apply the (binomial) crossover for DE algorithm
#'
#' Binomial crossover stipulates that crossover will occur on each 
#' of the D values in a solution whenever a randomly generated number
#' between 0 and 1 is within the CR range.
#'
#' @param x target ranking
#' @param v donor individual (mutaded x)
#' @param CR Crossover range
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo \email{giuliomazzeo@gmail.com}
#'
#' @return modified ranking
#' 
#' @importFrom stats runif
#' 
#' @export


crossover=function(x,v,CR){
  
  
  if (!is.matrix(x)){x=matrix(x,1,length(x))}
  
  D = ncol(x)
  
  u = matrix(0,1,D)
  
  for (i in 1:D){
    
    if (runif(1) > CR){
      
      u[i] = v[i]
    }
    else{
      
      u[i] = x[i]
      
    }
  }
  u
}