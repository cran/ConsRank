#' Stirling numbers of the second kind 
#'
#' Denote the number of ways to partition a set of n objects into k non-empty subsets
#'
#' @param n (integer): the number of the objects
#' @param k (integer <=n): the number of the non-empty subsets (buckets)
#'
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' S \tab  \tab the stirling number of the second kind\cr
#' SM \tab       \tab a matrix showing, for each k (on the columns) in how many ways the n objects (on the rows) can be partitioned}
#' 
#'
#' @references Comtet, L. (1974). Advanced Combinatorics: The art of finite and infinite expansions. D. Reidel, Dordrecth, The Netherlands.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' 
#' @examples
#' parts<-stirling2(4,2)
#'
#' @keywords Stirling numbers of second kind
#' 
#' @export




stirling2 <- function(n,k){
  #Stirling numbers of the second kind, denoting the
  #number of ways to partition a set of n objects into k non-empty subsets
  #
  #input:  n (integer): the number of the objects
  #        k (integer): the number of the non-empty subsets (buckets)
  #
  #output:
  #        S:  the stirling number of the second kind
  #        SM: the matrix showing, for each k (on the columns) in how many
  #            ways the n objects (on the rows) can be partitioned
  
  if (k==0){
    
    if (n==0){S<-1}else{S<-0}
    
    SM<-matrix(0,n,n)
    
  }  else  {
    
    SM<-matrix(NA,n,k)
    SM[,1]<-1
    
    ind <- (1:k-1)*n+1:k
    
    SM[ind]<-1
    
    for (i in 2:n){
      
      crit<-min((i-1),k)
      if (crit>=2){
        
        for (j in 2:crit){
          
          SM[i,j]<-SM[(i-1),(j-1)]+j*SM[(i-1),j]
          
        }
        
      }
      
    }
    
    S<-SM[n,k]
    
    
  } #end else
  
  
  return(list(S=S, SM=SM))
  
}

#---------------------------------------------------------
