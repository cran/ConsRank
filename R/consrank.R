#'Branch-and-bound and heuristic algorithms to find consensus (median) ranking according to the Kemeny's axiomatic approach
#'
#'Branch-and-bound, Quick , FAST and DECOR algorithms to find consensus (median) ranking according to the Kemeny's axiomatic approach. The median ranking(s) can be restricted to be necessarily a full ranking, namely without ties   
#'
#' @param X A n by m data matrix, in which there are n judges and m objects to be judged. Each row is a ranking of the objects which are represented by the columns. If X contains the rankings observed only once, the argument wk can be used
#' @param wk Optional: the frequency of each ranking in the data
#' @param ps If PS=TRUE, on the screen some information about how many branches are processed are displayed.
#' @param algorithm Specifies the used algorithm. One among "BB", "Quick", "FAST" and "DECOR". algorithm="BB" is the default option.
#' @param full Specifies if the median ranking must be searched in the universe of rankings including all the possible ties (full=FALSE) or in the restricted space of full rankings (permutations). full=FALSE is the default option.
#' @param itermax maximum number of iterations for FAST and DECOR algorithms. itermax=10 is the default option.
#' @param np For DECOR algorithm only: the number of population individuals. np=15 is the default option.
#' @param gl For DECOR algorithm only: generations limit, maximum number of consecutive generations without improvement. gl=100 is the default option.
#' @param ff For DECOR algorithm only: the scaling rate for mutation. Must be in [0,1]. ff=0.4 is the default option.
#' @param cr For DECOR algorithm only: the crossover range. Must be in [0,1]. cr=0.9 is the default option.
#' @param proc For BB algorithm only: proc=TRUE allows the branch and bound algorithm to work in difficult cases, i.e. when the number of objects is larger than 15 or 25. proc=FALSE is the default option
#' 
#' @details The BB algorithm can take long time to find the solutions if the number objects to be ranked is 
#' large with some missing (>15-20 if full=FALSE, <25-30 if full=TRUE). 
#' quick algorithm works with a large number of items to be ranked. The solution is quite accurate. 
#' fast algorithm works with a large number of items to be ranked by repeating several times the quick algorithm with different random starting points.
#' decor algorithm works with a very large number of items to be ranked. 
#' For decor algorithm, empirical evidence shows that the number of population individuals (the 'np' parameter) can be set equal to 10, 20 or 30
#' for problems till 20, 50 and 100 items. Both scaling rate and crossover ratio (parameters
#' 'ff' and 'cr') must be set by the user. The default options (ff=0.4, cr=0.9) work well
#' for a large variety of data sets
#' All algorithms allow the user to set the option 'full=TRUE' if the median ranking(s) must be searched in the restricted space of permutations instead of in the unconstrained universe of rankings of n items including all possible ties 
#' 
#'  
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}#' 
#' 
#' @examples
#' data(Idea)
#' RevIdea<-6-Idea 
#' # as 5 means "most associated", it is necessary compute the reverse ranking of 
#' # each rankings to have rank 1 = "most associated" and rank 5 = "least associated"
#' CR<-consrank(RevIdea)
#' CR<-consrank(RevIdea,algorithm="quick")
#' #CR<-consrank(RevIdea,algorithm="fast",itermax=10)
#' #not run
#' #data(EMD)
#' #CRemd<-consrank(EMD[,1:15],wk=EMD[,16],algorithm="decor",itermax=1)
#' #data(APAFULL)
#' #CRapa<-consrank(APAFULL,full=TRUE)
#' 
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' D'Ambrosio, A., Amodio, S., and Iorio, C. (2015). Two algorithms for finding optimal solutions of the Kemeny rank aggregation problem for full rankings. Electronic Journal of Applied Statistical Analysis, 8(2), 198-213.
#' Amodio, S., D'Ambrosio, A. and Siciliano, R. (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, 249(2), 667-676.
#' D'Ambrosio, A., Mazzeo, G., Iorio, C., and Siciliano, R. (2017). A differential evolution algorithm for finding the median ranking under the Kemeny axiomatic approach. Computers and Operations Research, vol. 82, pp. 126-138. 
#' 
#' 
#' 
#'
#' @keywords Consensus ranking
#' @keywords Median ranking
#' @keywords Branch-and-bound
#' @keywords Quick algorithm
#' @keywords Fast algorithm
#' @keywords Differential evolution
#' @keywords Genetic algorithms
#' 
#' @importFrom stats runif
#' @importFrom rlist list.rbind list.cbind
#' @importFrom methods is
#' @importFrom proxy dist
#' @importFrom gtools combinations
#' 
#' @export















consrank<-function(X,wk=NULL,ps=TRUE,algorithm="BB",full=FALSE,itermax=10,np=15,gl=100,ff=0.4,cr=0.9,proc=FALSE){
  # X:         a N by M data matrix in which there are N judges and M objects to be judged. 
  #            Each row is a ranking of the objects which are represented by the columns. 
  #            Alternatively X can contain the rankings observed only once in the sample. 
  #            In this case the argument Wk must be used
  # wk:	       Optional: the frequency of each ranking in the data
  # FULL:      Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space 
  #            of full rankings.
  # algorithm: One among "BB" (default, Emond and Mason's Branch and Bound algorithm), 
  #            "quick", "fast", (Quick algorithm by Amodio' D'Ambrosio and Siciliano)  and "decor"
  #            (Differential evolution algorithm by D'Ambrosio, Mazzeo, Iorio and Siciliano)
  # itermax:   maximum number of iteration (only if algorithm = "FAST" or algorithm = "Decor", defaul 10)
  # np:        Number of population individuals (only if algorithm="Decor", defaul 15)
  # gl:        Generations limit: maximum number of consecutive generations without improvement
  #            (only if algorithm="Decor", defaul 100)
  # ff:        The scaling rate for mutation. Must be in [0,1] (only if algorithm="Decor", defaul 0.4)
  # cr:        The crossover range. Must be in [0,1] (only if algorithm="Decor", defaul 0.9)
  
  
  if (is(X,"data.frame")) {
    #colnames(X)=NULL
    X<-as.matrix(X)
  }
  
  if (is(nrow(X),"NULL")){X<-matrix(X,nrow<-1)}
  
  if(algorithm=="BB"){ #if1 
    
    #-- start check 1: if there is a larger number of objects maybe it is better use another algorithm
    if (full==FALSE & proc==FALSE){
      
      if (nrow(X)>1 & ncol(X)>15){ 
        stop('the number of items is larger than 15. If data matrix contains several ties 
               or if there are few judges, choose another algorthm (Quick, FAST or DECOR). If you want use the Branch-and-bound algorithm, then 
               set proc="TRUE"')
        #end check 1
      } else {out <- EMConsn(X,Wk=wk,PS=ps)}#end condition full false   
      
      
    } else if (full==TRUE & proc==FALSE) { 
      
      #-- start check 2: if there is a larger number of objects maybe it is better use another algorithm
      
      if (nrow(X)>1 & ncol(X)>25){
        stop('the number of items is larger than 25. If data matrix contains several ties 
               or if there are few judges, choose another algorthm (Quick, FAST or DECOR). If you want use the Branch-and-bound algorithm, then 
               set proc="TRUE"')
        #end check 2
      } else { out <- BBFULLn(X,Wk=wk,PS=ps) }#end condition full true 
      
    } else if (full==FALSE & proc==TRUE){
      
      out <- EMConsn(X,Wk=wk,PS=ps)
      
    } else if (full==TRUE & proc==TRUE){
      
      out <- BBFULLn(X,Wk=wk,PS=ps)
      
    }
    
  } #end BB
    
    
    
 
  #-- Quick algorithm
  
  if (algorithm=="quick") {
    out <- QuickConsn(X,Wk=wk, FULL=full, PS=ps)
  } 
  
  #-- FAST algorithm
  
  if (algorithm=="fast"){
    out <- FASTconsn(X, Wk=wk, maxiter=itermax, FULL=full, PS=ps)
  } 
  
  #-- DECOR
  
  if (algorithm=="decor"){
    out<-FASTDECORn(X,Wk=wk,maxiter=itermax,NP=np,L=gl,FF=ff,CR=cr,FULL=full,PS=ps)
  }
  
  return(out)

}
#---------------------------------------------------------------------------------------------------------------------

breakties<-function(X){
  X<-rank2order(X)
  open<-"{"
  close<-"}"
  X<-gsub(open,"", X,fixed=TRUE)
  X<-gsub(close,"", X,fixed=TRUE)
  X<-order2rank(X)
  return(X)
}

#----------------------------------------------------------------------------------------------------------------------------

EMConsn <- function(X,Wk=NULL,PS=TRUE)  {
  #Emond and Mason Branch and Bound algorithm to find median ranking
  #X is a data matrix in which the rows are the judges and the columns indicates the objects
  #Wk is the vector of weigths
  
  if (is(X,"data.frame")) {
    #colnames(X)=NULL
    X<-as.matrix(X)
  }
  
  
  
  M <- nrow(X)
  N<-ncol(X)
  callps<-PS
  tic <- proc.time()[3]
  if (M==1) {
    consensus <- X
    TauX <- 1
  } else {
    if (!is(Wk,"NULL")) {
      
      if (is(Wk,"numeric")) {
        Wk<-matrix(Wk,ncol=1)
      }
      
      cij <- combinpmatr(X,Wk)
    } else {
      cij <- combinpmatr(X)
    }
    
    if (sum(cij==0)==length(cij)){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    if (sum(sign(cij+diag(N)))==length(cij)){
      print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
      return()
      
    }
    
    R<-findconsensusBB(cij)
    cons1<-BBconsensus(R,cij,FULL=FALSE,PS=FALSE)
    consensus1<-cons1$cons
    Po<-cons1$pen
    consensus<-BBconsensus2(consensus1,cij,Po,PS=callps,FULL=FALSE)
  }
  
  
  if (nrow(consensus)==1) {
    
    Sij<-scorematrix(consensus)
    
    if (!is(Wk,"NULL")){
      TauX<-sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      TauX<-sum(cij*Sij) / (  M*(N*(N-1)) )
    }
    
  } else {
    
    TauX<-matrix(0,nrow(consensus),1)
    
    for (k in 1:nrow(consensus)) {
      
      Sij<-scorematrix(t(matrix(consensus[k,])))
      
      if (!is(Wk,"NULL")) {
        
        TauX[k,1] <- sum(cij*Sij) / ( sum(Wk)*(N*(N-1)) )
        
      } else {
        
        TauX[k,1] <- sum(cij*Sij) / (M*(N*(N-1)))
        
      }
      
    }
    
  }
  toc <- proc.time()[3]
  colnames(consensus)<-colnames(X) 
  #consensus<-reordering(consensus)
  eltime<-toc-tic
  return(list(Consensus=reordering(consensus), Tau=TauX, Eltime=eltime) )
}

#--------------------------------------------------------------------------------------

findconsensusBB <- function(cij,FULL=FALSE) {
  
  X<-mat.or.vec(1,ncol(cij)) + 1
  N<-ncol(X)
  indici<-combinations(N,2)
  for (j in 1:nrow(indici)) {
    if ( sign(cij[indici[j,1],indici[j,2]]) == 1 & sign(cij[indici[j,2],indici[j,1]]) == -1 ) {
      X[indici[j,1]]<-X[indici[j,1]]+1
    } else if ( sign(cij[indici[j,1],indici[j,2]]) == -1 & sign(cij[indici[j,2],indici[j,1]]) == 1 ) {
      X[indici[j,2]]<-X[indici[j,2]] + 1
    } else if (sign(cij[indici[j,1],indici[j,2]]) == -1 & sign(cij[indici[j,2],indici[j,1]]) == -1 ) {
      X[indici[j,1]]<- NA
    } else if (sign(cij[indici[j,1],indici[j,2]]) == 1 & sign(cij[indici[j,2],indici[j,1]]) == 1 ){
      X[indici[j,1]]<-X[indici[j,1]]+1
      X[indici[j,2]]<-X[indici[j,2]] + 1
    }
    
  }
  
  X<-(N+1)-X;
  if (FULL==TRUE){
    X<-breakties(X)
  }
  return(X)
}

#--------------------------------------------------------------------------------------------------

BBconsensus <- function(RR,cij,FULL=FALSE,PS=FALSE) {
  
  #Branch and Bound Algorithm to find the the consensus ranking PART I As modified by D'AMBROSIO (2008).
  #Find the first approximation to the consensus ranking. Most of the time CR
  #is a solution, maybe not unique
  #Input:
  #       RR -> First solution candidate to be the consensus ranking
  #       cij -> Combined Input Matrix of the M individuals which judge n
  #       objects
  #Output:
  #       CR -> Consensus Ranking
  #
  #
  #References: Amodio et al.,2015; D'Ambrosio et al., 2016.
  
  CR<-RR
  sij<-scorematrix(RR)
  Po <- sum(abs(cij))-sum(cij*sij)
  a <- t(matrix(sort(RR,decreasing = TRUE)))
  ord <- t(matrix(order(RR,decreasing = TRUE)))
  R<-RR
  addpenalty<-matrix(0,length(a),1)
  
  # exploration of the initial solution
  for (k in 2:length(a)) {
    #print(k)
    b <- 1:k
    R <- ReorderingBB(R)
    KR<-t(matrix(R[ord[b]]))
    KR<-KR[-length(KR)]
    MO<-max(KR)
    MI<-min(KR)
    aa<-1
    KO<-1
    KR[length(KR)+1]<-MO+1
    R[ord[b]]<-KR
    candidate<-matrix(0,nrow(RR), ncol(RR))
    Pb <- matrix(0, 1, 1)
    while (KO==1)  {
      #browser()
      candidate<-rbind(candidate,R)
      #if (ncol(candidate>ncol(RR))) {
      
      #print(dim(candidate))
      #print(dim(R))
      #print(class(candidate))
      #print(class(R))
      
      # }
      if (aa==1){
        candidate<-matrix(candidate[-1,],1,ncol(candidate))
      }
      
      Sij<-scorematrix(matrix(candidate[aa,],1,ncol(R)))
      # print(Sij)
      #print(candidate)
      #flush.console()
      Pb<-rbind(Pb,sum(abs(cij))-sum(cij*Sij))
      if (aa==1) {
        Pb<-matrix(Pb[-1,],1,1)
      }
      # print(Pb)
      if (Pb[aa]==0) {
        
        CR <- R
        Po <- 0
        Pc <- 0
        break
      }
      Pc<-1
      if(FULL==TRUE){
        R[ord[b[length(b)]]] <- R[ord[b[length(b)]]]-2 }else{
          R[ord[b[length(b)]]] <- R[ord[b[length(b)]]]-1
        }
      if (MI-R[ord[b[length(b)]]] > 1) {
        KO <- 0
      }
      aa<-aa+1
      
    }
    
    if (PS==TRUE) {
      
      dsp2<-paste("evaluated",nrow(candidate),"branches",sep=" ")
      print(dsp2)
      
    }
    
    if (Pc == 0) {
      break
    }
    
    minp<-min(Pb)
    posp<-which(Pb==min(Pb))
    
    if (minp<=Po) {
      Po<-minp
      CR<-t(matrix(candidate[posp[1],]))
      R<-CR
      addpenalty[k,1]<-PenaltyBB2(cij,R,ord[b])
    } else {
      R <- CR
      addpenalty[k,1]<-PenaltyBB2(cij,R,ord[b])
    }
    
    candidate <- mat.or.vec(nrow(R), ncol(R))
    Pb <- mat.or.vec(1, 1)
    
  }
  
  if (Pc==0) {
    Po<-0
    addpenalty <- 0
  }  else {
    Poo<-sum(addpenalty)
  }
  
  SIJ <- scorematrix(CR)
  Po<-sum(addpenalty)
  
  return(list(cons=CR,pen=Po))
}

#------------------------------------------------------------------------------------------------------------------------------

BBconsensus2 <- function(RR,cij,Po,PS=TRUE,FULL=FALSE) {
  ##Core code for the computation of the consensus ranking. Branch-and-bound 
  ##algorithm by Emond and Mason
  ForFULL<-FULL
  CR<-RR
  a <- t(matrix(sort(RR)))
  ord <- t(matrix(order(RR)))
  r<-ReorderingBB(RR)
  BR.R<-r #initialize ranking
  BR.P<-0 #initialize penalty
  WCi<-1
  lambda<-1
  nobj<-ncol(RR)
  while (WCi == 1){
    if (PS==TRUE) {
      dsp1<-paste("round",lambda,sep=" ")
      print(dsp1) 
    }
    
    for (k in 2:ncol(a)) { #primary loop: add the k^th better object ranked
      
      B <- nrow(BR.R)    #branches size
      #      print(B)
      #      flush.console()
      
      b<-1:k
      
      
      for (nb in 1:B) { #Secondary loop:  check the branches created by "nb"
        
        BR.R[nb,] <- ReorderingBB(t(matrix(BR.R[nb,])))
        rpbr<-branches(matrix(BR.R[nb,],1,),cij,b,Po,ord,matrix(BR.P[nb]),FULL=ForFULL)
        R<-rpbr$cR
        Pbr<-rpbr$pcR     
        #        print(nrow(R))
        #        flush.console()
        
        
        if (is.null(R)) {
          
          #if (nb==1) {
          #
          #    JR=0
          #    
          #    next
          #    
          #} else {     
          
          next
          
          #}
        } else {         #process good rankings
          
          if (nb==1) {     #starting point
            
            #JR=nrow(R)
            
            KR.R<-R
            KR.P<-Pbr
            
            
            
          } else {   #it is not the starting point
            
            KR.R<-rbind(KR.R,R)
            KR.P<-rbind(KR.P,Pbr)
            
            
          }
          
        }
        
        #JR = nrow(KR.R)  #update size of branches
        
      } #end secondary loop
      
      if (is.null(R)) {
        
        if (nb==B & nb!=1) { #If at the last iteration of the secondary loop all of thye rankings are not good 
          
          rm(BR.R) #BR.R = NULL    #rm(BR.R)
          rm(BR.P) #BR.P = NULL    #rm(BR.P)
          BR.R <- KR.R
          BR.P <- KR.P
          KR.R <- NULL
          KR.P <- NULL
          #rm(JR) #JR = NULL
          
        } else {
          
          next
          
        }
        
      } else {
        
        rm(BR.R) #BR.R = NULL  #rm(BR.R)
        rm(BR.P) #BR.P = NULL  #rm(BR.P)
        BR.R <- KR.R
        BR.P <- matrix(KR.P)
        KR.R <- NULL
        KR.P <- NULL        
        #        KR.R = NULL
        #        KR.P = NULL
        #rm(JR)  #JR = NULL
        
      }
      
      if (PS==TRUE) {
        
        dsp2<-paste("evaluating",B,"branches",sep=" ")
        print(dsp2)
        
      }
      
    } #end primary loop
    
    #AccPen = BR.P
    
    
    SSP<-matrix(which(BR.P==min(BR.P)))
    MinP<-min(BR.P)
    PenMin<-Po-MinP
    
    if (PenMin==0) {
      
      #CR=t(matrix(BR.R[SSP,]))
      CR<-matrix(BR.R[SSP,],length(SSP),nobj)
      WCi <- 0
      
    } else {
      
      Po<-MinP
      WCi<-1
      lambda<-lambda+1
      #nRR<-t(matrix((BR.R[SSP[1],])))
      nRR<-matrix((BR.R[SSP[1],]),1,nobj)
      rm(BR.R)
      rm(BR.P)
      BR.R<-nRR
      BR.P<-0
      a <- t(matrix(sort(BR.R)))
      ord <- t(matrix(order(BR.R)))
      rm(nRR)
      
    }
    
  }   #end while
  
  CR
}

#--------------------------------------------------------------------------------------

branches <- function(brR,cij,b,Po,ord,Pb,FULL=FALSE) {
  
  candidate <- findbranches(brR,ord,b,FULL)
  Pb <- matrix( rep(Pb,nrow(candidate)))
  
  CR<-mat.or.vec(nrow(candidate),ncol(candidate))
  addpenalty<-matrix(0,nrow(candidate),1)
  QR<-mat.or.vec(nrow(candidate),ncol(candidate))
  
  for (gm in 1:nrow(candidate)) {
    
    CR[gm,]<-candidate[gm,]
    addpenalty[gm,]<-PenaltyBB2(cij,candidate[gm,],ord[b])
    
    if ( (Pb[gm]+addpenalty[gm,]) > Po) {
      
      #      CR[gm,]=-10.0e+15
      #      addpenalty[gm]=-10.0e+15
      CR[gm,]<-NA
      addpenalty[gm,]<-NA      
      
    }
    QR[gm,]<-CR[gm,]
  }
  Pbr<-addpenalty+Pb
  #  idp=Pbr<0
  idp<-which(is.na(Pbr))
  
  if (length(idp)==0) {
    
    R<-QR
    
  } else if (length(idp)==nrow(QR)) {
    
    Pbr<-NULL
    Pb<-NULL
    R<-NULL
    
  } else {
    #    Pbr=t(matrix(Pbr[idp==FALSE,],1))
    #    if (sum(idp==F)==1) {
    #      R=t(matrix(QR[idp==FALSE,]))
    #    } else {
    #      R=QR[idp==FALSE,]
    #    }
    
    Pbr<-matrix(Pbr[-idp],length(Pbr[-idp]),1)
    R<-QR[-idp,]
    if (is(nrow(R),"NULL")){R<-matrix(R,1,length(R))}
    
  }
  
  return(list(cR=R,pcR=Pbr))
}

#---------------------------------------------------------------------------------------

Penalty <- function(CR,cij,indice)   #indice must be order(CR)
{
  if (CR[indice[1,1]] < CR[indice[1,2]]) { #case 1, the first object is preferred
    #to the second object
    if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
      Po<-0
    } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po <- cij[indice[1,2],indice[1,1]]-cij[indice[1,1],indice[1,2]]  #     cji-cij
    } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po <- cij[indice[1,2],indice[1,1]] #cji
    }
  } else if (CR[indice[1,1]] > CR[indice[1,2]]) { #case 2 the first object is not
    #preferred to the second one
    if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
      Po<- cij[indice[1,1],indice[1,2]]-cij[indice[1,2],indice[1,1]]  #cij-cji
    } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po <- 0
    } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po <- cij[indice[1,1],indice[1,2]] #cij
    }
  } else if (CR[indice[1,1]] == CR[indice[1,2]]) { #case 3 they are in a tie
    if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
      Po <- -cij[indice[1,2],indice[1,1]] #-cj
    } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po <- -cij[indice[1,1],indice[1,2]] #-cij
    } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1 |
               sign(cij[indice[1,1],indice[1,2]]) == 0 & sign (cij[indice[1,2],indice[1,1]]) == 0 )  {
      Po <- 0
    }
    
  }
  Po
}

#-----------------------------------------------------------------------------------------------------------

PenaltyBB2 <- function(cij,candidate,ord)   #indice must be order(CR)
  
  ## DETERMINATION OF PENALTIES FOR THE BRANCH AND BOUND ALGORITHM
  
{
  
  Ds<-t(mat.or.vec(1,(length(ord)-1)));
  addpenalty<-t(mat.or.vec(1,(length(ord)-1)));
  
  for (k in 1:(length(ord)-1)) {
    
    Ds[k,1]<-sign(candidate[ord[length(ord)]]-candidate[ord[k]]);
    
    if (Ds[k,1]==1) {
      
      
      if ( sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
        addpenalty[k,1]<-cij[ord[length(ord)],ord[k]]-cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 & sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 & sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 & sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]<-cij[ord[length(ord)],ord[k]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 & sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]<-0
      }
    }
    else if (Ds[k,1]==-1) {
      if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
        addpenalty[k,1]<-0
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]<-cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]<-cij[ord[k],ord[length(ord)]]-cij[ord[length(ord)],ord[k]]
      }
    }
    
    else if (Ds[k,1]==0) {
      if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]])  == -1 ) {
        addpenalty[k,1] <- -cij[ord[k],ord[length(ord)]]
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 1  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 1 && sign(cij[ord[k],ord[length(ord)]]) == 0  ||
                  sign(cij[ord[length(ord)],ord[k]]) == 0 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]<-0
      } else if ( sign(cij[ord[length(ord)],ord[k]]) == -1 && sign(cij[ord[k],ord[length(ord)]]) == 1) {
        addpenalty[k,1]=-cij[ord[length(ord)],ord[k]]
      }
    }
  }
  
  addpenalty<-sum(addpenalty)
  
}

#---------------------------------------------------------------------------------------------------------------------

ReorderingBB <- function(RR) {
  
  RR <- RR+1
  R <- RR;
  k <- ncol(R)
  neword <- order(R)
  indexing <- mat.or.vec(1, ncol(R)-1)
  for (j in (k-1):1) {
    indexing[j] <- R[neword[j+1]]-R[neword[j]]
  }
  
  if (sum(indexing==0)>0) {
    J <- 1
    while (J<=ncol(indexing)) {
      if (indexing[J]==0) {
        R[neword[J+1]]<-R[neword[J]]
        J<-J+1
      }
      else if (indexing[J]>0) {
        R[neword[J+1]] <- R[neword[J]]+2
        J<-J+1
      }
    }
  }
  else  {
    J <- 1
    while (J<= ncol(indexing)) {
      R[neword[J+1]] <- R[neword[J]] + 2
      J<-J+1
    }
  }
  R
}

#-------------------------------------------------------------------------------


findbranches <- function(R,ord,b,FULL=FALSE) {
  
  KR<-t(matrix(R[ord[b]]))
  KR<-KR[-length(KR)]
  MO<-max(KR)
  MI<-min(KR)
  aa<-1
  KO<-1
  KR[length(KR)+1]<-MO+1;
  R[ord[b]]<-KR
  candidate<-mat.or.vec(nrow(R), ncol(R))
  
  while (KO==1)  {
    candidate<-rbind(candidate,R)
    
    if (aa==1){
      candidate<-matrix(candidate[-1,],1,ncol(candidate))
    }
    if (FULL==FALSE){
      R[ord[b[length(b)]]]<-R[ord[b[length(b)]]]-1 }else{
        R[ord[b[length(b)]]]<-R[ord[b[length(b)]]]-2
      }
    
    if (MI-R[ord[b[length(b)]]] > 1) {
      
      KO<-0
      
    }
    
    aa<-aa+1
    
  }
  
  Rt<-candidate
  
}

#---------------------------------------------------------------------------------------------------

QuickConsn <- function(X,Wk=NULL, FULL=FALSE,PS=FALSE)   {
  
  
  if (is(X,"data.frame")) {
    #colnames(X)=NULL
    X<-as.matrix(X)
  }
  
  callfull<-FULL
  callps<-PS
  
  M <- nrow(X)
  N<-ncol(X)
  tic <- proc.time()[3]
  
  if (M==1) {
    consensus <- X
    TauX <- 1
  } else {
    if (!is(Wk,"NULL")) {
      
      if (is(Wk,"numeric")) {
        Wk<-matrix(Wk,ncol=1)
      }
      
      cij <- combinpmatr(X,Wk)
    } else {
      cij <- combinpmatr(X)
    }
    
    if (sum(cij==0)==length(cij)){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    if (sum(sign(cij+diag(N)))==length(cij)){
      print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
      return()
      
    }
    
    R<-findconsensusBB(cij,FULL=callfull)
    R1<-(N+1)-R
    consensusA <- BBconsensus(R,cij,FULL=callfull,PS=callps)$cons
    consensusB <- BBconsensus(consensusA,cij,FULL=callfull,PS=callps)$cons
    consensusC <- BBconsensus(R1,cij,FULL=callfull,PS=callps)$cons
    consensusD <- BBconsensus(consensusC,cij,FULL=callfull,PS=callps)$cons
    consensus <- unique(reordering(rbind(consensusA,consensusB,consensusC,consensusD)))
    howcons <- nrow(consensus)
    
    
  }
  #d=kemenyd(X,consensus$cons)
  
  Taux<-matrix(0,nrow(consensus),1)
  for (k in 1:nrow(consensus)) {
    
    #Sij=scorematrix(t(as.matrix(consensus[k,])))
    Sij<-scorematrix(matrix(consensus[k,],1,N))
    
    if (!is(Wk,"NULL")){
      Taux[k,1]<-sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      Taux[k,1]<-sum(cij*Sij) / (  M*(N*(N-1)) )
    }
    
  }
  
  if (howcons>1) {
    nco<-which(Taux==max(Taux))
    if (length(nco)>1) {
      consensus<-consensus[nco,]
      Taux<-matrix(rep(max(Taux),nrow(consensus)),nrow(consensus),1)
    } else {
      Taux<-max(Taux)
      #consensus <- t(matrix(consensus[nco,]))
      consensus <- matrix(consensus[nco,],1,N)
    }
  }
  colnames(consensus)<-colnames(X) 
  toc <- proc.time()[3]
  eltime<-toc-tic
  return(list(Consensus=reordering(consensus), Tau=Taux, Eltime=eltime) )
}

#------------------------------------------------------------------------------------

FASTconsn <- function(X, Wk=NULL, maxiter=50, FULL=FALSE, PS=FALSE)   {
  
  
  if (is(X,"data.frame")) {
    #colnames(X)=NULL
    X<-as.matrix(X)
  }
  
  M <- nrow(X)
  N<-ncol(X)
  
  tic <- proc.time()[3]
  if (M==1) {
    CR <- X
    Taux <- 1
  } else {
    if (!is(Wk,"NULL")) {
      
      if (is(Wk,"numeric")) {
        Wk<-matrix(Wk,ncol=1)
      }
      
      cij <- combinpmatr(X,Wk)
    } else {
      cij <- combinpmatr(X)
    }
    
    if (sum(cij==0)==length(cij)){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    if (sum(sign(cij+diag(N)))==length(cij)){
      print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
      return()
      
    } 
    
    CR<-matrix(0,maxiter,ncol(X))
    for (iter in 1:maxiter) {
      
      if (FULL==TRUE){
        R<-matrix(sample(1:ncol(X)),1,ncol(X))
      } else {
        
        if (iter%%2==0) {
          R<-matrix(sample(1:ncol(X),replace=TRUE),1,ncol(X))
        } else {
          R<-matrix(sample(1:ncol(X)),1,ncol(X))
        }
      }
      
      consensus1 <- BBconsensus(R,cij, FULL)
      cons<-matrix(consensus1$cons,1,ncol(X))
      consensus <- BBconsensus(cons,cij, FULL)
      #print(cons)
      #print(R)
      #flush.console()
      CR[iter,]<-matrix(consensus$cons,1,ncol(X))
      if (PS==TRUE) {
        
        dsp1<-paste("Iteration",iter,sep=" ")
        print(dsp1)
      }
      
    }
  }
  #d=kemenyd(X,consensus$cons)
  
  Taux<-matrix(0,nrow(CR),1)
  for (k in 1:nrow(CR)) {
    Sij<-scorematrix(matrix(CR[k,],1,ncol(X)))
    if (!is.null(Wk)){
      Taux[k,]<-sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      Taux[k,]<-sum(cij*Sij) / (  M*(N*(N-1)) )
    }
  }
  
  CR<-reordering(CR)
  indice<-which(Taux==max(Taux));
  Taux<-max(Taux)
  CR<-matrix(CR[indice,],ncol=N)
  if (nrow(CR>1)){
    CR<-unique(CR)
  }
  if (!is(dim(CR),"NULL")) {
    Taux<-matrix(rep(Taux,nrow(CR)))
  }
  
  colnames(CR)<-colnames(X) 
  toc <- proc.time()[3]
  eltime<-toc-tic
  return(list(Consensus=CR, Tau=Taux, Eltime=eltime) )
}

#------------------------------------------------------------------------------

DECOR <- function(X,Wk=NULL,NP=15,L=100,FF=0.4,CR=0.9,FULL=FALSE){
  
  #check if X is a matrix
  if (is(X,"data.frame")) {
    #colnames(X)=NULL
    X<-as.matrix(X)
  }
  
  M <- nrow(X)
  N<-ncol(X)
  tic <- proc.time()[3]  
  
  #check if there are trivial solutions
  
  if (M==1) { 
    consensus <- X
    TauX <- 1
    
  } else {
    
    if (!is(Wk,"NULL")) {
      
      if (is(Wk,"numeric")) {
        Wk<-matrix(Wk,ncol=1)
        NJ<-sum(Wk)
      }
      
      cij <- combinpmatr(X,Wk)
      
    } else {
      
      cij <- combinpmatr(X)
      NJ<-nrow(X)
    }
    
    
    if (sum(cij==0)==length(cij)){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    if (sum(sign(cij+diag(N)))==length(cij)){
      print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
      return()
      
    }
    
    
    COR<-DECORcore(cij,NJ,NP,L,FF,CR,FULL)
    
    Consensus<-COR$ConsR
    TauX<-COR$Tau
    colnames(Consensus)<-colnames(X) 
    row.names(Consensus)<-NULL
    
    
    
    
  }
  toc<-proc.time()[3]
  eltime<-toc-tic
  return(list(Consensus=reordering(Consensus), Tau=TauX, Eltime=eltime) )
}

#-------------------------------------------------------------------------------------------------

DECORcore <- function(cij,NJ,NP=15,L=50,FF=0.4,CR=0.9,FULL=FALSE){
  
  # DECoR Differential Evolution for COnsensus Ranking
  
  # Discrete version of Differential Evolution algorithm
  # specifically created for solving the Consensus Ranking problem,
  # AKA social choice, AKA rank aggregation.
  # This file uses external functions for mutation and crossover during
  # the genetic evolutions of population.
  # Another function is needed to discretize and correct child solutions
  # when mutation and crossover generate out-of-bound values and duplicates.
  #
  # Input parameters:
  #
  # NP            - The number of population individuals
  #
  # L             - Generations limit: maximum number of consecutive 
  #                 generations without improvement
  #
  # F             - The scaling rate for mutation. Must be in [0,1]
  #
  # CR            - The crossover range. Must be in [0,1]
  #
  # cij           - Combined input matrix (see Emond and Mason B&B)
  #
  # NJ            - Number of judges (needed for cost computation)
  #
  # FULL          - FULL = 1 search in the space of full rankings
  #
  # Output parameters:
  #
  # besti         - Best Individuals:
  #                 matrix of best individuals for every generation
  #
  # bestc         - Best Costs:
  #                 vector of best individuals' cost for every gen
  #
  # bests         - Best Solutions:
  #                 matrix with "all" the best solutions (founded)
  #                 for checking if more than one optimal solution is found
  #
  # Notes:
  #
  # mutation      - rand/1/bin > best/1/bin
  # 
  #Authors : Giulio Mazzeo and Antonio D'Ambrosio
  
  
  # preparation
  
  tic <- proc.time()[3]
  N<-nrow(cij)        # number of objects                     
  costs       <- matrix(0,1,NP)  # array of initial costs
  
  # initialize the population (random selection)
  population <- matrix(0,(NP-1),N)
  
  for (k in 1:(NP-1)){ population[k,] <- sample(N)}
  
  
  # insert a very good candidate
  
  population<-rbind(population,findconsensusBB(cij,FULL))
  
  #  if (FULL==TRUE){
  #    population[NP,]=order(population[NP,])
  #  }
  
  
  
  # compute costs of initial population
  costs<-matrix(0,NP,1)
  taos<-costs
  for (i in 1:NP){
    
    COTA<-combincost(population[i,],cij,NJ)
    costs[i]<-COTA$cp
    taos[i]<-COTA$tp
  }
  
  # store the best individual and cost of initial population
  bestc <- min(costs)
  bestind <- which(costs==min(costs))
  bestT <- max(taos)
  besti<-population[bestind,]
  
  
  # generation index
  g <- 2
  
  # evolve for generations
  no_gain <- 0
  
  while (no_gain < L){
    
    
    # individuals mutation
    for (i in 1:NP){
      
      
      # apply mutation
      evolution <- mutaterand1(population,FF,i);
      
      # apply crossover
      evolution <- crossover(population[i,],evolution,CR)
      
      
      # apply discretization and convert to rank
      
      if (FULL==TRUE){
        
        evolution<-order(evolution)}
      
      else{
        
        evolution <- childtie(evolution)
      }
      
      # apply selection, hold the best individual
      COTAN <- combincost(evolution,cij,NJ)
      cost_new<-COTAN$cp
      ta_new<-COTAN$tp
      
      
      if (cost_new < costs[i]){
        population[i,] <- evolution
        costs[i] <- cost_new
        taos[i]<-ta_new
      }
      
    }
    
    # store the best individual of current generation
    
    bestco <- min(costs)
    bestc<-rbind(bestc,bestco)
    bestind <- which.min(costs)
    bestTa <- max(taos)
    bestT<-rbind(bestT,bestTa)
    bestin<-population[bestind,]
    besti<-rbind(besti,bestin)
    
    
    # check if this generation improved solutions
    if (bestc[g] == bestc[(g-1)]){
      
      no_gain <- no_gain + 1}
    
    else{
      
      no_gain <- 0
    }
    
    
    # next generation
    g <- g + 1
    
  } #end while
  
  # select ALL the best solutions
  
  indexes <- which(bestc==min(bestc))
  if (FULL==TRUE){ #if1
    
    if (length(indexes)==1){ #if2
      bests<-childclosint(matrix(besti[indexes,],1,N))}
    else{
      bests<-matrix(0,length(indexes),N)
      for (j in 1:length(indexes)){
        bests[j,]<-childclosint(besti[indexes[j],])
      } #end for
    } #end if2
    
  } else { #if FULL = FALSE
    
    if(length(indexes)==1){
      
      bests <- reordering(matrix(besti[indexes,],1,N))
      
    } else {
      
      bests <- reordering(besti[indexes,])}
    
  } #end if1
  
  avgTau <- bestT[indexes]
  
  ConsR<-unique(bests)
  Tau<-matrix(rep(avgTau,nrow(ConsR)),nrow(ConsR),1)
  
  
  toc <- proc.time()[3]
  eltime<-toc-tic
  return(list(ConsR=ConsR,Tau=Tau,besti=besti,bestc=bestc,bests=bests,avgTau=avgTau,bestT=bestT,Eltime=eltime))
  
}

#----------------------------------------------------------------------------------------------------------------

FASTDECORn <- function(X,Wk=NULL,maxiter=10,NP=15,L=100,FF=0.4,CR=0.9,FULL=FALSE,PS=TRUE){
  
  #check if X is a matrix
  if (is(X,"data.frame")) {
    #colnames(X)=NULL
    X<-as.matrix(X)
  }
  
  M <- nrow(X)
  N<-ncol(X)
  tic <- proc.time()[3]  
  
  #check if there are trivial solutions
  
  if (M==1) { 
    consensus <- X
    TauX <- 1
    
  } else {
    
    if (!is(Wk,"NULL")) {
      
      if (is.numeric(Wk)) {
        Wk<-matrix(Wk,ncol=1)
        NJ<-sum(Wk)
      }
      
      cij <- combinpmatr(X,Wk)
      
    } else {
      
      cij <- combinpmatr(X)
      NJ<-nrow(X)
    }
    
    
    if (sum(cij==0)==length(cij)){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    if (sum(sign(cij+diag(N)))==length(cij)){
      print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
      return()
      
    }
    
    
    sol<-matrix(0,1,N)
    taos<-0
    
    for (iter in 1:maxiter){
      COR<-DECORcore(cij,NJ,NP,L,FF,CR,FULL)
      
      sol<-rbind(sol,COR$ConsR)
      taos<-rbind(taos,COR$Tau)
      
      if (PS==TRUE) {
        
        if (iter%%5==0){
          
          dsp1<-paste("Iteration",iter,sep=" ")
          print(dsp1)
        }
        
      }
      
      
    }
    
  }
  
  sol<-sol[-1,]
  taos<-matrix(taos[-1],(length(taos)-1),1)
  
  if (is(nrow(sol),"NULL")){sol<-matrix(sol,1,N)}
  
  bestindex<-which(taos==max(taos))
  sol<-sol[bestindex,]
  taos<-max(taos)
  
  if (is(nrow(sol),"NULL")){
    sol<-matrix(sol,1,N)
    Consensus<-sol
  } else {
    Consensus<-unique(sol)
  }
  
  if (is(nrow(Consensus),"NULL")){
    Consensus<-matrix(Consensus,1,N)
  }
  
  TauX<-matrix(rep(taos),nrow(Consensus),1)
  colnames(Consensus)<-colnames(X) 
  row.names(Consensus)<-NULL
  
  
  toc<-proc.time()[3]
  eltime<-toc-tic
  return(list(Consensus=reordering(Consensus), Tau=TauX, Eltime=eltime) )
}

#----------------------------------------------------------------------


childclosint<-function(r){
  
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
  
  
  D <- length(r)
  
  # closest integer
  x <- round(r)
  
  # correct out of bound
  for (i in 1:D){
    
    
    if (x[i]>D | x[i]<1){
      
      r <- sample(D)
      
      x[i] <- r[1]
      
    }
    
    
    # correct duplicates
    
    
    
    C <- setdiff(union(x,1:D),intersect(x,1:D))
    
    if (length(C)==0){
      
      x<-x
      
    }
    else
    {
      U <- unique(x)
      id <- which(!duplicated(x)) 
      ix <- setdiff(union(id,1:D),intersect(id,1:D))
      x[ix] <- C[sample(length(ix))]
      
    } #end else
    
  } #end for
  x
} #end function

#----------------------------------------------------------------

childtie <- function( r ){
  
  #CHILDTIE Summary of this function goes here
  #   Detailed explanation goes here
  o <- matrix(rank(r,ties.method="average"),1,ncol(r))
  o <- reordering(o)
  o
  
}


#--------------------------------------------------------------

crossover<-function(x,v,CR){
  
  
  if (!is(x,"matrix")){x<-matrix(x,1,length(x))}
  
  D <- ncol(x)
  
  u <- matrix(0,1,D)
  
  for (i in 1:D){
    
    if (runif(1) > CR){
      
      u[i] <- v[i]
    }
    else{
      
      u[i] <- x[i]
      
    }
  }
  u
}

#---------------------------------------------------------------

mutaterand1<-function(X,FF,i){
  
  
  
  D <- nrow(X)
  
  a <- sample(D)
  
  for (j in 1:3){
    
    if (a[j]==i){
      
      a[j]<-a[4]
      
    }
    
  }
  
  r1 <- a[1]
  r2 <- a[2]
  r3 <- a[3]
  
  # apply mutation (Rand1Bin)
  v <- X[r1,] + FF*(X[r2,]-X[r3,])
  v
}


#----------------------------------------------------------------



BBFULLn <- function(X,Wk=NULL,PS=TRUE)  {
  #Branch and Bound algorithm to find median ranking in the space of full rankings
  #X is a data matrix in which the rows are the judges and the columns indicates the objects
  #Wk is the vector of weigths
  if (is(X,"data.frame")) {
    #colnames(X)=NULL
    X<-as.matrix(X)
  }
  
  
  
  M <- nrow(X)
  N<-ncol(X)
  tic <- proc.time()[3]
  if (M==1) {
    consensus <- X
    TauX <- 1
  } else {
    if (!is.null(Wk)) {
      
      if (is(Wk,"numeric")) {
        Wk<-matrix(Wk,ncol=1)
      }
      
      cij <- combinpmatr(X,Wk)
    } else {
      cij <- combinpmatr(X)
    }
    
    if (sum(cij==0)==nrow(cij)^2){
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking")
      return()
      
    } 
    
    R<-findconsensusBB(cij,FULL=TRUE)
    cons1<-BBconsensus(R,cij,FULL=TRUE)
    consensus1<-cons1$cons
    Po<-cons1$pen
    consensus<-BBconsensus2(consensus1,cij,Po,PS,FULL=TRUE)
  }
  
  
  if (nrow(consensus)==1) {
    
    Sij<-scorematrix(consensus)
    
    if (!is(Wk,"NULL")){
      TauX<-sum(cij*Sij) / ( sum(Wk)* (N*(N-1)) )
    } else {
      TauX<-sum(cij*Sij) / (  M*(N*(N-1)) )
    }
    
  } else {
    
    TauX<-matrix(0,nrow(consensus),1)
    
    for (k in 1:nrow(consensus)) {
      
      Sij<-scorematrix(t(matrix(consensus[k,])))
      
      if (!is(Wk,"NULL")) {
        
        TauX[k,1] <- sum(cij*Sij) / ( sum(Wk)*(N*(N-1)) )
        
      } else {
        
        TauX[k,1] <- sum(cij*Sij) / (M*(N*(N-1)))
        
      }
      
    }
    
  }
  toc <- proc.time()[3]
  colnames(consensus)<-colnames(X) 
  #consensus<-reordering(consensus)
  eltime<-toc-tic
  return(list(Consensus=reordering(consensus), Tau=TauX, Eltime=eltime) )
}

#-------------------------------------------------------------------------------

combincost <- function(ranking,cij,M){
  
  
  
  if (!is(ranking,"matrix")){ranking<-matrix(ranking,1,length(ranking))}
  
  N <- ncol(ranking)
  sij <- scorematrix(ranking)
  
  # max distance
  maxdist <- (N*(N-1))
  
  # computing the tau
  t <- sum(sum(cij*sij))/(M*(N*(N-1)))
  
  # computing the distance (from tau)
  c <- (maxdist/2)*(1-t)*M;    
  
  return(list(tp=t,cp=c))
}






  

  