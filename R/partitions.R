#'Generate partitions of n items constrained into k non empty subsets
#'
#'Generate all possible partitions of n items constrained into k non empty subsets. It does not generate the universe of rankings constrained into k buckets. 
#'
#' @param n a (integer) number denoting the number of items
#' @param k The number of the non-empty subsets. Default value is NULL, in this case all the possible partitions are displayed
#' @param items items: the items to be placed into the ordering matrix. Default are the first c small letters
#' @param itemtype to be used only if items is not set. The default value is "L", namely letters. Any other symbol produces items as the first c integers
#' 
#' @details If the objects to be ranked is large (>15-20) with some missing, it can take long time to find the solutions. If the searching space is 
#' limited to the space of full rankings (also incomplete rankings, but without ties), use the function BBFULL or the functions FASTcons and QuickCons 
#' with the option FULL=TRUE.
#' 
#' @return the ordering matrix (or vector)
#' 
#' @examples
#' X<-partitions(4,3)
#' #shows all the ways to partition 4 items (say "a", "b", "c" and "d" into 3 non-empty subets
#'  #(i.e., into 3 buckets). The Stirling number of the second kind (4,3) indicates that there
#'  #are 6 ways.
#' s2<-stirling2(4,3)$S
#' X2<-order2rank(X) #it transform the ordering into ranking
#' 
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @seealso \code{\link{stirling2}} Stirling number of second kind.
#' @seealso \code{\link{rank2order}} Convert rankings into orderings. 
#' @seealso \code{\link{order2rank}} Convert orderings into ranks.
#' @seealso \code{\link{univranks}} Generate the universe of rankings given the input partition
#'
#' 
#' @importFrom rlist list.rbind list.cbind
#' 
#' @export





partitions <- function(n,k=NULL,items=NULL,itemtype="L"){
  
  #generate all possible partitions of n items constrained into k non empty subsets
  #input: 
  #n:         a (integer) number denoting the number of items
  #k:         the number of the non-empty subsets. Default value is NULL,
  #           in this case all the possible partitions are displayed
  #items: the items to be placed into the ordering matrix. Default are the
  #       first c small letters
  #itemtype:to be used only if items is not set. The default value is "L", namely
  #         letters. Any other symbol produces items as the first c integers
  #
  #OUTPUT: the ordering matrix (or vector)
  #
  #example
  # X=partitions(4,3) #shows all the ways to partition 4 itemns (say "a", "b", "c" and "d" into 3 non-empty subets
  #                   #(i.e., into 3 buckets). The Stirling number of the second kind (4,3) indicates that there
  #                   #are 6 ways. The output is
  #                   #> X
  #                   #     [,1] [,2] [,3] [,4]
  #                   #[1,] "{a" "b}" "c"  "d" 
  #                   #[2,] "{a" "c}" "b"  "d" 
  #                   #[3,] "a"  "{b" "c}" "d" 
  #                   #[4,] "{a" "d}" "b"  "c" 
  #                   #[5,] "a"  "{b" "d}" "c" 
  #                   #[6,] "a"  "b"  "{c" "d}"
  #
  # X2=order2rank(X) #it transform the ordering into ranking   
  #                  #> X2
  #                  #     a b c d
  #                  #[1,] 1 1 2 3
  #                  #[2,] 1 2 1 3
  #                  #[3,] 1 2 2 3
  #                  #[4,] 1 2 3 1
  #                  #[5,] 1 2 3 2
  #                  #[6,] 1 2 3 3
  
  
  #check
  
  if(!is(k,"NULL")){
    if(k>n) stop('k must be less or equal to n')
    if(k<=0 || (k%%1)!=0) stop('k must be positive integers')
  }
  
  if(n<=0 || n%%1 !=0) stop('n must be a positive integer')
  #end checks
  
  if(is(items,"NULL")){
    
    if(itemtype=="L"){
      
      items<-letters[seq(1:n)]
      
    } else {
      
      items<-as.character(seq(1:n))
    }
  }
  
  if(is(k,"NULL")){
    parts<-partall2(n,items)
  } else {
    parts<-partk(n,k,items)
  }
  
  
  return(parts)  
  
}

#-----------------------------------------------------

insert <- function(part,n,elements=NULL,stack=1){
  
  # perform all possible insertions
  
  #----preprocessing----------------------------------
  if(is(elements,"NULL")){elements<-letters[seq(1:n)]}
  
  if (length(part)==1){ #if1
    
    l<-1 
    open<-NULL
    close<-NULL
    
  } else {
    
    if (!is(part,"matrix")){part<-matrix(part,1,length(part))}
    
    
    open<-unlist(gregexpr(pattern ='\\{',part))
    close<-unlist(gregexpr(pattern ='\\}',part))
    
    if (length(which(open>=1)) + length(which(close>=1)) == 0){ #if2
      
      
      l <- ncol(part)
      
    } else {
      
      sequ<-rowSums(cbind(apply(part,1,function(x) unlist(gregexpr(pattern ='\\{',part))),
                         apply(part,1,function(x) unlist(gregexpr(pattern ='\\}',part)))))
      
      #sequ=rbind(sequ,seq(1:length(sequ)))
      
      buckets<-vector(mode="list")
      
      checka<-FALSE
      
      i<-1
      k<-1
      inda<-NA
      while(checka==FALSE){
        
        
        if (sequ[i]==0){  #if the sequence starts with a tied ranking
          checkb<-FALSE
          while(checkb==FALSE){
            
            if (sequ[i]==1){
              
              inda<-cbind(inda,i)
              inda<-inda[-1]
              buckets[[k]]<-part[inda]
              inda<-NA
              
              k<-k+1
              i<-i+1
              checkb<-TRUE
              
            } else {
              
              
              inda<-cbind(inda,i)
              i<-i+1
              
              
            }
            
            
          } #end while checkb
          
        } else if (sequ[i]==-2){ # if the sequence starts with an object not in a tie
          checkc<-FALSE
          while(checkc==FALSE){
            inda<-cbind(inda,i)
            inda<-inda[-1]
            buckets[[k]]<-part[inda]
            i<-i+1
            k<-k+1
            checkc<-TRUE
          }
          
          
          
          # if (i==length(sequ)){
          #   inda=cbind(inda,i)
          #   inda=inda[-1]
          #   buckets[[k]]=part[inda]
          #   i=i+1
          #   checkc=TRUE
          # } else if (sequ[(i+1)]==0){
          #   
          #   inda=cbind(inda,i)
          #   inda=inda[-1]
          #   buckets[[k]]=part[inda]
          #   inda=NA
          #   i=i+1
          #   k=k+1
          #   checkc=TRUE
          #   
          # } else {
          #   
          #   
          #   inda=cbind(inda,i)
          #   i=i+1
          #   
          # } 
          
          # } #end while c
          
          
          
          
        } # end if statements
        
        
        if (i>length(sequ)){checka<-TRUE}
        
        
      } #while 1
      
      l<-length(buckets)
      
    }  #end if2
    
    
  }#end if1
  
  
  
  #---------end preprocessing---------------------------
  
  en <- elements[n]
  if (stack==0){
    m<-l
  } else if (stack==2){
    m <- 1
  } else {
    m<-l+1
  }
  
  
  part_i <- vector(mode = "list", length = m)
  part_ib <- vector(mode = "list", length = m)
  
  
  if (length(which(open>=1)) + length(which(close>=1)) == 0){ #if there are no ties
    
    if(stack<=1){#if stack <=1
      
      #for (j in 1:(m-1)){part_i[[j]]=part}
      
      for (i in 1:l){
        
        part_i[[i]]<-part
        
        if (is(part_i[[i]],"numeric")){part_i[[i]]<-as.character(part_i[[i]])}
        
        indice<-1:length(part_i[[i]])
        indice[indice>=i]<-indice[indice>=i]+1
        part_ib[[i]][indice]<-part_i[[i]]
        part_ib[[i]][i]<-paste("{",part_i[[i]][i],sep="")
        part_ib[[i]][(i+1)]<-paste(en,"}",sep="")
      }
      
    }#end if stack <=1
    
    if (stack>=1){
      
      part_ib[[m]]<-as.character(cbind(part,en))
      
    }
    
  } else  {  #if there are ties
    
    if (stack<=1){ #if stack >=1
      
      for (i in 1:l){
        
        if (is(buckets[[i]],"numeric")){buckets[[i]]<-as.character(buckets[[i]])}
        
        buck<-buckets[[i]]
        buck<-gsub("\\{", "", buck)
        buck<-gsub("\\}", "", buck)
        buck<-cbind(list.cbind(buck),en)
        buck[1]<-paste("{",buck[1],sep="")
        buck[length(buck)]<-paste(buck[length(buck)],"}",sep="")
        rasp<-buckets
        rasp[[i]]<-buck
        part_ib[[i]]<-unlist(rasp)
        
      }
      
    }#end if stack >=1
    
    if (stack>=1){
      part_ib[[m]]<-as.character(cbind(list.cbind(unlist(buckets)),en))
    }
    
    
  } #end if there are no ties
  
  
  return(list.rbind(part_ib))
  
}



#-----------------------------------

partall2 <- function(n, elements=NULL){
  
  
  # Return the cell list of all partitions of the integer set {1:n}
  # Output LIST is the cell of size (b x 1), where b is Bell number Bn
  # Each list{j} is a partition of {1:n}
  
  if(is(elements,"NULL")){elements<-letters[seq(1:n)]}
  
  if (n==1){
    
    lista<-elements[n]
    
  } else {
    
    lp<-partall2((n-1),elements) #recursive call to partall2
    
    if(is(nrow(lp),"NULL")){k<-length(lp)}else{k<-nrow(lp)}
    
    for (i in 1:k){
      
      if(k==1){llp<-lp[i]}else{llp<-lp[i,]}
      
      pa<-insert(llp,n,elements)
      if (i==1){lista<-pa}else{lista<-rbind(lista,pa)}
      
      
    }
    
  }
  
  return(lista)
}

#------------------------------------

partk<-function(n,k,elements=NULL){
  
  if(is(elements,"NULL")){elements<-letters[seq(1:n)]}
  
  
  if(k==1){
    lista<-matrix(elements,1,length(elements))
    lista[1]<-paste("{",lista[1],sep="")
    lista[length(lista)]<-paste(lista[length(lista)],"}",sep="")
  } else if (k==length(elements)) {
    lista<-matrix(as.character(elements),1,length(elements))
  } else {
    
    
    
    m<-n-k+1
    
    L<-vector(mode="list", length=m)
    
    for (j in 1:m){
      
      L[[j]]<-elements[1:j]
      
      if(j>1){
        L[[j]][1]<-paste("{",L[[j]][1],sep="")
        L[[j]][length(L[[j]])]<-paste(L[[j]][length(L[[j]])],"}",sep="")
      }
      
    }
    
    #SN=stirling2(n,k)
    
    for (kk in 2:k){ #for #2
      
      L[[1]]<-elements[1:kk]
      
      for (j in 2:m){# for #3
        
        nu<-j+kk-1
        
        #lista=vector(mode="list", length= SN$SM[nu,kk])
        
        lp<-L[[j]]
        
        if(is(nrow(lp),"NULL")){
          tt<-1} else {
            tt<-nrow(lp)
          }
        
        if(tt>1){
          for(i in 1:tt){
            part_i<-insert(lp[i,],nu,elements,stack=2)
            if (i==1){lista1<-part_i}else{lista1<-rbind(lista1,part_i)}
          }
        } else {
          
          lista1<-insert(matrix(lp,1,length(lp)),nu,elements,stack=2)
        }
        
        lp<-L[[(j-1)]]
        
        if(is(nrow(lp),"NULL")){
          tt<-1} else {
            tt<-nrow(lp)
          }
        
        
        if(tt>1){
          
          for(i in 1:tt){
            part_i<-insert(lp[i,],nu,elements,stack=0)
            if (i==1){lista2<-part_i}else{lista2<-rbind(lista2,part_i)}
          }
        } else {
          lista2<-insert(matrix(lp,1,length(lp)),nu,elements,stack=0)
        }
        
        L[[j]]<-rbind(lista1,lista2)
        
      }#end for3
      
    }#end for2
    
    
    
    lista<-L[[m]]
  }#end first if
  
  return(lista)
  
  
}#end function

#------------------------------------------

