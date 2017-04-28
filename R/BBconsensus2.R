#' Core function in computing consensus ranking as defined by Emond and Mason (2002)
#'
#' Core function in computing consensus ranking as defined by Emond and Mason (2002), recalled by EMCons function
#'
#' @param RR A ranking
#' @param cij combined input matrix
#' @param Po current penalty
#' @param PS If PS=true, it prints the evaluating branches
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. FULL=TRUE if the function is called by BBFULL algorithm.
#'
#' @return median ranking
#' 
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' #
#' D'Ambrosio, A., Amodio, S., and Iorio, C. (2015). Two algorithms for finding optimal solutions of the Kemeny rank aggregation problem for full rankings. Electronic Journal of Applied Statistical Analysis, 8(2), 198-213.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @seealso \code{\link{EMCons}} Emond and Mason branch-and-bound algorithm
#' @seealso \code{\link{BBFULL}} D'Ambrosio et al. branch-and-bound algorithm for full rankings
#' 
#' @export

BBconsensus2 = function(RR,cij,Po,PS=TRUE,FULL=FALSE) {
  ##Core code for the computation of the consensus ranking. Branch-and-bound 
  ##algorithm by Emond and Mason
  ForFULL=FULL
  CR=RR
  a = t(matrix(sort(RR)))
  ord = t(matrix(order(RR)))
  r=ReorderingBB(RR)
  BR.R=r #initialize ranking
  BR.P=0 #initialize penalty
  WCi=1
  lambda=1
  nobj=ncol(RR)
  while (WCi == 1){
    if (PS==TRUE) {
      dsp1=paste("round",lambda,sep=" ")
      print(dsp1) 
    }
    
    for (k in 2:ncol(a)) { #primary loop: add the k^th better object ranked
      
      B = nrow(BR.R)    #branches size
#      print(B)
#      flush.console()
      
      b=1:k
      
      
      for (nb in 1:B) { #Secondary loop:  check the branches created by "nb"
        
        BR.R[nb,] = ReorderingBB(t(matrix(BR.R[nb,])))
        rpbr=branches(matrix(BR.R[nb,],1,),cij,b,Po,ord,matrix(BR.P[nb]),FULL=ForFULL)
        R=rpbr$cR
        Pbr=rpbr$pcR     
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
            
            KR.R=R
            KR.P=Pbr
            
            
            
          } else {   #it is not the starting point
            
            KR.R=rbind(KR.R,R)
            KR.P=rbind(KR.P,Pbr)
            
            
          }
          
        }
        
        #JR = nrow(KR.R)  #update size of branches
        
      } #end secondary loop
      
      if (is.null(R)) {
        
        if (nb==B & nb!=1) { #If at the last iteration of the secondary loop all of thye rankings are not good 
          
          rm(BR.R) #BR.R = NULL    #rm(BR.R)
          rm(BR.P) #BR.P = NULL    #rm(BR.P)
          BR.R = KR.R
          BR.P = KR.P
          KR.R = NULL
          KR.P = NULL
          #rm(JR) #JR = NULL
          
        } else {
          
          next
          
        }
        
      } else {
        
        rm(BR.R) #BR.R = NULL  #rm(BR.R)
        rm(BR.P) #BR.P = NULL  #rm(BR.P)
        BR.R = KR.R
        BR.P = matrix(KR.P)
        KR.R = NULL
        KR.P = NULL        
#        KR.R = NULL
#        KR.P = NULL
        #rm(JR)  #JR = NULL
        
      }
      
      if (PS==TRUE) {
        
        dsp2=paste("evaluating",B,"branches",sep=" ")
        print(dsp2)
        
      }
      
    } #end primary loop
    
    #AccPen = BR.P
    
    
    SSP=matrix(which(BR.P==min(BR.P)))
    MinP=min(BR.P)
    PenMin=Po-MinP
    
    if (PenMin==0) {
      
      #CR=t(matrix(BR.R[SSP,]))
      CR=matrix(BR.R[SSP,],length(SSP),nobj)
      WCi = 0
      
    } else {
      
      Po=MinP
      WCi=1
      lambda=lambda+1
      #nRR=t(matrix((BR.R[SSP[1],])))
      nRR=matrix((BR.R[SSP[1],]),1,nobj)
      rm(BR.R)
      rm(BR.P)
      BR.R=nRR
      BR.P=0
      a = t(matrix(sort(BR.R)))
      ord = t(matrix(order(BR.R)))
      rm(nRR)
      
    }
    
  }   #end while
  
  CR
}
