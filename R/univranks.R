#'Generate the universe of rankings
#'
#'Generate the universe of rankings given the input partition
#'
#' @param X A ranking, an ordering, a matrix of rankings, a matrix of orderings or a number
#' @param k Optional: the number of the non-empty subsets. It has to be used only if X is anumber. The default value is NULL, In this case the universe of rankings with n=X items are computed
#' @param ordering The universe of rankings must be returned as orderings (default) or rankings?
#' 
#' @details The function should be used with small numbers because it can generate a large number of permutations. The use of X greater than 9, of X matrices with more than 9 columns as input is not reccomended. 
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Runiv \tab  \tab The universe of rankings\cr
#' Cuniv \tab       \tab A list containing:\cr
#'  \tab R   \tab The universe of rankings in terms of rankings;\cr
#'  \tab Parts \tab for each ranking in input the produced rankings\cr
#'  \tab Univinbuckets \tab the universe of rankings within each bucket}
#' 
#' @examples
#' S2<-stirling2(4,4)$SM[4,] #indicates in how many ways 4 objects
#'                          #can be placed, respectively, into 1, 2,
#'                          #3 or 4 non-empty subsets.
#' CardConstr<-factorial(c(1,2,3,4))*S2  #the cardinality of rankings 
#'                                      #constrained into 1, 2, 3 and 4
#'                                      #buckets
#' Card<-sum(CardConstr)  #Cardinality of the universe of rankings with 4
#'                       #objects                                                              
#' U<-univranks(4)$Runiv #the universe of rankings with four objects
#'                      # we know that the universe counts 75 
#'                      #different rankings
#' Uk<-univranks(4,2)$Runiv    #the universe of rankings of four objects 
#'                            #constrained into k=2 buckets, we know they are 14
#' Up<-univranks(c(1,4,3,1))$Runiv  #the universe of rankings with 4 objects
#'                                 #for which the first and the fourth item
#'                                 #are tied
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @seealso \code{\link{stirling2}} Stirling number of second kind.
#' @seealso \code{\link{rank2order}} Convert rankings into orderings. 
#' @seealso \code{\link{order2rank}} Convert orderings into ranks.
#' @seealso \code{\link{partitions}} Generate partitions of n items constrained into k non empty subsets.
#'
#' 
#' @export



univranks = function(X,k=NULL,ordering=TRUE){
  if (!is(X,"numeric")){ 
    X<-order2rank(X) } else if (is(X,"numeric")){
      if (length(X)==1){X<-order2rank(partitions(X,k))}
    }
  univ<-allranksnew(X)
  if (ordering==TRUE){
    univ2<-rank2order(univ$R)
  } else {
    univ2<-univ$R
  }
  return(list(Runiv=univ2,Cuniv=univ))
}



#--------------------------------------





allranksnew <- function (X){
  
  #generate all possible rankings of n items given
  #the ranking X
  # allranks(order2rank(partall2(4)))
  
  
  if(is(nrow(X),"NULL")){
    r<-1
    c<-length(X)
    X<-matrix(X,r,c)
    nbuckts<-sortnbuckts<-length(unique(X)) #number of buckets in the ranking X
    ordnbuckts<-1
    #    walk=1
  } else {
    r<-nrow(X)
    c<-ncol(X)
    nbuckts<-apply(X,1,function(x) length(unique(x))) #buckets for each ranking in X
    sortnbuckts<-sort(nbuckts)
    ordnbuckts<-order(nbuckts)
    X<-X[ordnbuckts,]
    #    walk=c(0,diff(sortnbuckts))
    #    sortnbuckts=c(NA,sort(nbuckts))
  }
  
  
  
  
  ref<-seq(1:c)
  out<-matrix(NA,1,c)
  SS<-NULL
  allperms<-vector(mode="list")
  permsinbuck<-vector(mode="list")
  nb<-rowSums(X)
  
  
  buckid<-1
  for (i in 1:r){
    
    if(nb[i]==c){
      allbuck<-X[i,]
      SS<-rbind(SS,-1)
    } else if(nb[i]==(c*(c+1)/2)) {
      allbuck<-permutations(c,c)
      SS<-rbind(SS,-1)
    } else {
      
      ord <- rank(X[i,])
      orders <- tapply(ref, ord, sort)
      SS<-rbind(SS,length(orders))
      if(i==1){
        indici<-permutations(length(orders),length(orders))
      } else {
        if(SS[i]!=SS[(i-1)]){
          indici<-permutations(length(orders),length(orders))
          #SS=rbind(SS,length(orders))
        }
      }
      allbuck<-matrix(NA,nrow(indici),c)
      for (j in 1:nrow(indici)){
        for (k in 1:length(orders)){
          id<-unlist(orders[indici[j,k]])
          allbuck[j,id]<-k
        }
      }
    } #end else
    
    out<-rbind(out,allbuck)
    allperms[[i]]<-allbuck
    
  }
  
  out<-out[-1,]
  row.names(out)<-NULL
  
  
  
  css<-c(0,cumsum(table(sortnbuckts)))
  
  for (h in 1:(length(css)-1)){
    permsinbuck[[h]]<-list.rbind(allperms[((css[h]+1):css[(h+1)])])
  }
  
  return(list(R=out,Parts=allperms,Univinbucket=permsinbuck))
}