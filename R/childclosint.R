#' Transform the vector into ranking for DECOR
#'
#' Closest integer approach is used: the elements of x are rounded
#' to the closest integer. Then check if any solution exists outside 
#' of the bounds (and get it back inside the bounds randomly).
#' Finally repair the solution if repetitions exist.
#' 
#' @param r mutated vector
#'
#' @return a valid ranking
#'
#' @references Davendra, D., and Onwubolu, G. (2007). Enhanced differential evolution hybrid scatter search for discrete optimization. In Evolutionary Computation, 2007. CEC 2007. IEEE Congress on (pp. 1156-1162). IEEE.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo \email{giuliomazzeo@gmail.com}
#' 
#' @export


childclosint=function(r){
  
  # CHILDCLOSINT Transform the vector into ranking for DECoR
  #
  #  Closest integer approach is used: the elements of x are rounded
  # to the closest integer. Then check if any solution exists outside 
  # of the bounds (and get it back inside the bounds randomly).
  # Finally repair the solution if repetitions exist.
  #
  # Random approach is proven to be effective in:
  #   "Enhanced Differential Evolution..." by Davendra and Onwubolu
  #
  # $Author: Giulio Mazzeo $    $Email: giuliomazzeo@gmail.com $
  
  
  D = length(r)
  
  # closest integer
  x = round(r)
  
  # correct out of bound
  for (i in 1:D){
    
    
    if (x[i]>D | x[i]<1){
      
      r = sample(D)
      
      x[i] = r[1]
      
    }
    
    
    # correct duplicates
    
    
    
    C = setdiff(union(x,1:D),intersect(x,1:D))
    
    if (length(C)==0){
      
      x=x
      
    }
    else
    {
      U = unique(x)
      id = which(!duplicated(x)) 
      ix = setdiff(union(id,1:D),intersect(id,1:D))
      x[ix] = C[sample(length(ix))]
      
    } #end else
    
  } #end for
  x
} #end function
