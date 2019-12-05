#' Given a rank, it is transformed to a ordering 
#'
#' From ranking to ordering. IMPORTANT: check which symbol denotes tied rankings in the X matrix 
#' 
#'
#' @param X A ordering or a matrix containing orderings
#' @param TO symbol indicating the start of a set of items ranked in a tie
#' @param TC symbol indicating the end of a set of items ranked in a tie
#' @param items items to be placed into the ordering matrix. Default are the
#       first c small letters
#' @param itemtype to be used only if items=NULL. The default value is "L", namely
#         letters. Any other symbol produces items as the first c integers
#' 
#' @return a ordering or a matrix of orderings:
#' \tabular{lll}{
#' out \tab  \tab ranking or matrix of rankings}
#' 
#' @examples
#' data(APAred)
#' ord<-rank2order(APAred)
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export



rank2order <- function (X,items=NULL,TO="{",TC="}",itemtype="L"){
  
  #Given a ranking, it is transformed into an ordering 
  #input: 
  #X:     a n by c ranking matrix or a 1 by c rank vector
  #TO:    symbol that denotes the beginning of a block of ties to be placed. 
  #       Any symbol can be used, default is "{".
  #TC:    symbol that denotes the end of a block of ties to be insert. 
  #       Any symbol can be used, default is "}".  
  #items: the items to be placed into the ordering matrix. Default are the
  #       first c small letters
  #itemtype:to be used only if items is not set. The default value is "L", namely
  #         letters. Any other symbol produces items as the first c integers
  #
  #OUTPUT: the ordering matrix (or vector)
  
  
  if(is(nrow(X),"NULL")){
    r<-1
    c<-length(X)
    X<-matrix(X,r,c)
  } else {
    r<-nrow(X)
    c<-ncol(X)
  }
  
  if(is(items,"NULL")){
    
    if(itemtype=="L"){
      
      items<-letters[seq(1:c)]
      
    } else {
      
      items<-as.character(seq(1:c))
    }
  }
  
  out<-matrix(0,r,c)
  
  for (i in 1:r){
    
    ord <- rank(X[i,])
    orders <- tapply(items, ord, sort)
    check<-FALSE
    j<-1
    h<-1
    
    while(check==FALSE) {
      
      if (length(orders[[h]]) > 1) {
        
        k<-length(orders[[h]])
        nams<-matrix(orders[[h]],1,k)
        nams[1]<-paste(TO, nams[1], sep="")
        nams[k]<-paste(nams[k], TC, sep="")
        passo<-seq(j,(j+k-1))
        out[i,passo]<-nams
        j<-passo[k]+1
        h<-h+1
      } else {
        out[i,j] <- paste(orders[[h]])#, collapse = " ")
        j<-j+1
        h<-h+1
      }
      
      if(h>length(orders)){check<-TRUE}
      
    } #end for j
    
  } #end for i
  
  return(out)
}