#' Given an ordering, it is transformed to a ranking 
#'
#' From ordering to rank. IMPORTANT: check which symbol denotes tied rankings in the X matrix 
#' 
#'
#' @param X A ordering or a matrix containing orderings
#' @param TO symbol indicating the start of a set of items ranked in a tie
#' @param TC symbol indicating the end of a set of items ranked in a tie
#'
#' 
#' @return a ranking or a matrix of rankings:
#' \tabular{lll}{
#' R \tab  \tab ranking or matrix of rankings}
#' 
#' @examples
#' data(APAred)
#' ord=rank2order(APAred) #transform rankings into orderings
#' ran=order2rank(ord) #transform the orderings into rankings
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @export






order2rank=function(X,TO="{",TC="}"){
  #Given an ordering, it is transformed to a ranking 
  #input: 
  #X:  an ordering matrix or an ordering vector
  #TO: symbol that denotes the beginning of a block of ties. 
  #    Any symbol can be used, default is "{".
  #TC: symbol that denotes the end of a block of ties. 
  #    Any symbol can be used, default is "}".  
  #
  #IMPORTANT: check which symbol denotes tied rankings
  #
  #OUTPUT: the ranking matrix (or vector)
  
  Xl<-X #duplicate the ordering
  
  if (is(nrow(X),"NULL")){
    r<-1
    c<-length(X)
  } else {
    r<-nrow(X)
    c<-ncol(X)
  }
  
  
  # openb=unlist(gregexpr(pattern ='\\{',Xl))
  # closeb=unlist(gregexpr(pattern ='\\}',Xl))
  Xl<-gsub(TO, "", Xl,fixed=TRUE)
  Xl<-gsub(TC, "", Xl,fixed=TRUE)
  
  if (r==1){items<-as.character(sort(Xl))}else{items<-as.character(sort(Xl[1,]))}
  R<-matrix(NA,r,c)
  colnames(R)<-items
  Rref<-seq(1:c)
  
  for (i in 1:r){
    if (r==1){ #if x is an ordering vector
      openb<-unlist(gregexpr(pattern = TO, X, fixed=TRUE))
      closeb<-unlist(gregexpr(pattern = TC, X, fixed=TRUE))
      
      if(sum(rowSums(cbind(openb,closeb))==-2)==c){ #if there are no ties
        
        for(j in 1:c){
          
          R[which(items==X[j])]<-j
          
        }
        
      } else { #if there are ties
        check<-FALSE
        pos<-1
        iter<-1
        id<-rowSums(cbind(openb,closeb))
        ido<-which(id==0)
        idc<-which(id==1)
        iterid<-1
        stp<-0
        
        j<-1
        while(check==FALSE){
          R[which(items==Xl[j])]<-pos
          
          if (stp==1){ido<-rep(0,iterid)}
          
          if (j==ido[iterid]){
            
            for (i in ido[iterid]:idc[iterid]){
              
              R[which(items==Xl[i])]<-pos
              #if (i==idc[iterid]){pos=pos+1}
              
            }
            
            pos<-pos+1
            j<-idc[iterid]+1
            iterid<-iterid+1
            
          } else {
            
            pos<-pos+1
            j<-j+1
            
          }
          
          
          #if(( j>c) || (iterid>length(idc) )){check=TRUE}
          if( j>c){check<-TRUE}
          if(iterid>length(idc)){stp<-1}
          
          
        } #end while
        
        
        
      }
      
    } else { #if X is an ordering matrix
      
      openb<-unlist(gregexpr(pattern =TO, X[i,], fixed=TRUE))
      closeb<-unlist(gregexpr(pattern =TC, X[i,], fixed=TRUE))
      
      if(sum(rowSums(cbind(openb,closeb))==-2)==c){ #if there are no ties
        
        for(j in 1:c){
          
          R[i,which(items==X[i,j])]<-j
          
        }
        
      } else { #if there are ties
        check<-FALSE
        pos<-1
        iter<-1
        id<-rowSums(cbind(openb,closeb))
        ido<-which(id==0)
        idc<-which(id==1)
        iterid<-1
        stp<-0
        
        j<-1
        while(check==FALSE){
          R[i,which(items==Xl[i,j])]<-pos
          
          if (stp==1){ido<-rep(0,iterid)}
          
          if (j==ido[iterid]){
            
            for (k in ido[iterid]:idc[iterid]){
              
              R[i,which(items==Xl[i,k])]<-pos
              #if (i==idc[iterid]){pos=pos+1}
              
            }
            
            pos<-pos+1
            j<-idc[iterid]+1
            iterid<-iterid+1
            
          } else {
            
            pos<-pos+1
            j<-j+1
            
          }
          
          
          #if(( j>c) || (iterid>length(idc) )){check=TRUE}
          if( j>c){check<-TRUE}
          if(iterid>length(idc)){stp<-1}
          
          
        } #end while
        
        
        
      }
      
    }# end if X is a matrix
    
    
  } #end principal loop (for i=1:r)
  
  return(R)
  
}#end function

#----------------------------------------------------------------------------------