#' Auxiliary function called by DECORcore
#'
#' Calculates the sum of distances between a candidate ranking and the data set of rankings
#'
#' @param ranking the candidate ranking
#' @param cij combined input matrix
#' @param M number of judges
#'
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' tp \tab  \tab tauX associated to the ranking\cr
#' cp \tab       \tab cost associated to the ranking}
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo
#' 
#' @export


combincost = function(ranking,cij,M){
  

  
  if (!is.matrix(ranking)){ranking=matrix(ranking,1,length(ranking))}
  
  N = ncol(ranking)
  sij = scorematrix(ranking)
  
  # max distance
  maxdist = (N*(N-1))
  
  # computing the tau
  t = sum(sum(cij*sij))/(M*(N*(N-1)))
  
  # computing the distance (from tau)
  c = (maxdist/2)*(1-t)*M;    
  
  return(list(tp=t,cp=c))
}
