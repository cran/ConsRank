#' Auxiliary function
#'
#' Apply discretization and convert to rank
#'
#' @param r a candidate rank vector
#' 
#' @return a ranking
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo \email{giuliomazzeo@gmail.com}
#' 
#' @export

childtie = function( r ){
  
  #CHILDTIE Summary of this function goes here
  #   Detailed explanation goes here
  o = matrix(rank(r,ties.method="average"),1,ncol(r))
  o = reordering(o)
  o
  
}