#' Mutation phase
#'
#' Creates mutation vector v based on the current population X. We use the
#' rand/1/bin system
#'
#' @param X population matrix
#' @param FF scaling factor
#' @param i population index to be ignored
#' 
#' @return the mutated vector
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo \email{giuliomazzeo@gmail.com}
#' 
#' @export

mutaterand1=function(X,FF,i){
  
  
  
  D = nrow(X)
  
  a = sample(D)
  
  for (j in 1:3){
    
    if (a[j]==i){
      
      a[j]=a[4]
      
    }
    
  }
  
  r1 = a[1]
  r2 = a[2]
  r3 = a[3]
  
  # apply mutation (Rand1Bin)
  v = X[r1,] + FF*(X[r2,]-X[r3,])
  v
}
