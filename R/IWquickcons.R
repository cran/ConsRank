## Item weighted QuickCons

#' The item-weighted Quick algorithm to find up to 4 solutions to the consensus ranking problem
#'
#' The item-weighted Quick algorithm finds up to 4 solutions. Solutions reached are most of the time optimal solutions.
#'
#' @param X A N by M data matrix in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once in the sample. In this case the argument Wk must be used
#' @param w A M-dimensional row vector (individually weighted items), or a M by M matrix (item similarities)
#' @param Wk Optional: the frequency of each ranking in the data
#' @param full Default full=FALSE. If full=TRUE, the searching is limited to the space of full rankings.
#' @param PS Default PS=FALSE. If PS=TRUE the number of evaluated branches is diplayed
#'
#' @details The item-weigthed Quick algorithm finds up the consensus (median) ranking according to the Kemeny's axiomatic approach. The median ranking(s) can be restricted to be necessarily a full ranking, namely without ties.
#'
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' Consensus \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged item-weighted TauX rank correlation coefficient\cr
#' Eltime\tab   \tab Elapsed time in seconds}
#'
#' @examples
#' #Individually weighted items
#'data("German")
#'w=c(10,5,5,10)
#'iwquickcons(X= German,w= w)
#'
#' #Item similirity weights
#'data(sports)
#'dim(sports)
#'P=matrix(NA,nrow=7,ncol=7)
#'P[1,]=c(0,5,5,10,10,10,10)
#'P[2,]=c(5,0,5,10,10,10,10)
#'P[3,]=c(5,5,0,10,10,10,10)
#'P[4,]=c(10,10,10,0,5,5,5)
#'P[5,]=c(10,10,10,5,0,5,5)
#'P[6,]=c(10,10,10,5,5,0,5)
#'P[7,]=c(10,10,10,5,5,5,0)
#'iwquickcons(X= sports, w= P)
#'
#' @author Alessandro Albano \email{alessandro.albano@unipa.it} \cr
#' Antonella Plaia \email{antonella.plaia@unipa.it}
#' @references Amodio, S., D'Ambrosio, A. and Siciliano, R. (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, 249(2), 667-676. \cr
#' Albano, A. and Plaia, A. (2021).  Element weighted Kemeny distance for  ranking data. Electronic  Journal  of  Applied Statistical Analysis, doi: 10.1285/i20705948v14n1p117
#'
#' @seealso \code{\link{consrank}}
#'
#' @keywords Item-weighted Quick algorithm
#'
#' @export
#'
iwquickcons<-function (X,w, Wk = NULL, full = FALSE, PS = FALSE)
{
  if(!is(dim(w),"NULL")){
    
    X <- t(apply(X,1,rank,ties.method ="min"))
    ws_k<- w
    
  } else{
    
    X <- X[,which(w!=0)]
    X <- t(apply(X,1,rank,ties.method ="min"))
    w <- w[which(w!=0)]
    # original ws_k <- matrix(,length(w),length(w))
    ws_k <- matrix(NA,length(w),length(w))
    
    for (j in 1:length(w)) {
      for (i in 1:length(w)){
        ws_k[j,i] <- w[j]*w[i]
      }}
  }
  
  if (class(X)[1] == "data.frame") {
    X = as.matrix(X)
  }
  callfull = full
  callps = PS
  M = nrow(X)
  N = ncol(X)
  tic = proc.time()[3]
  if (M == 1) {
    consensus = X
    TauX = 1
  }
  else {
    if (!is.null(Wk)) {
      if (is.numeric(Wk)) {
        Wk = matrix(Wk, ncol = 1)
      }
      cij = iwcombinpmatr(X,w, Wk)
    }
    else {
      cij = iwcombinpmatr(X,w)
    }
    #trivial solution 1
    if (sum(cij == 0) == length(cij)) {
      print("Combined Input Matrix contains only zeros: any ranking in the reference universe is a median ranking. The output contains the identity permutation")
      cr <- matrix(seq(1:N),nrow=1)
      colnames(cr) <- colnames(X)
      return(list(Consensus=cr, Tau=0, Eltime=NA))
    }
    #trivial solution 2
    if (full==FALSE){
      
      if (sum(sign(cij + diag(N))) == length(cij)) {
        print("Combined Input Matrix contains only positive values: the median ranking is the all-tie solution")
        cr <- matrix(rep(1,N),nrow=1)
        Sij = scorematrix(cr)
        if (!is(Wk,"NULL")) {
          Taux <- sum(cij * Sij)/(sum(Wk)*(sum(ws_k)-sum(diag(ws_k))))
        }
        else {
          Taux <- sum(cij * Sij)/(M*(sum(ws_k)-sum(diag(ws_k))))
        }
        colnames(cr) <- colnames(X)
        return(list(Consensus=cr, Tau=Taux, Eltime=NA))
      }
    }
    
    R = findconsensusBB(cij)
    R1 = (N + 1) - R
    consensusA = BBconsensus(R, cij, FULL=callfull,PS=callps)$cons
    consensusB = BBconsensus(consensusA, cij, FULL=callfull,PS=callps)$cons
    consensusC = BBconsensus(R1, cij, FULL=callfull,PS=callps)$cons
    consensusD = BBconsensus(consensusC, cij, FULL=callfull,PS=callps)$cons
    consensus = unique(reordering(rbind(consensusA, consensusB,
                                        consensusC, consensusD)))
    howcons = nrow(consensus)
  }
  
  
  Taux = matrix(0, nrow(consensus), 1)
  for (k in 1:nrow(consensus)) {
    Sij = scorematrix(matrix(consensus[k, ], 1, N))
    if (!is(Wk,"NULL")) {
      Taux[k, 1] = sum(cij * Sij)/(sum(Wk)*(sum(ws_k)-sum(diag(ws_k))))
    }
    else {
      Taux[k, 1] = sum(cij * Sij)/(M*(sum(ws_k)-sum(diag(ws_k))))
    }
  }
  if (howcons > 1) {
    nco = which(Taux == max(Taux))
    if (length(nco) > 1) {
      consensus = consensus[nco, ]
      Taux = matrix(rep(max(Taux), nrow(consensus)), nrow(consensus),
                    1)
    }
    else {
      Taux = max(Taux)
      consensus = matrix(consensus[nco, ], 1, N)
    }
  }
  colnames(consensus) = colnames(X)
  toc = proc.time()[3]
  eltime = toc - tic
  return(list(Consensus = reordering(consensus), Tau = Taux,
              Eltime = eltime))
}

#-------------------


#' Item-weighted Combined input matrix of a data set
#' 
#' Compute the item-weighted Combined input matrix of a data set as defined by Albano and Plaia (2021)
#' 
#' @param X A data matrix N by M, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
#' @param w A M-dimensional row vector (individually weighted items), or a M by M matrix (item similarities)
#' @param Wk Optional: the frequency of each ranking in the data
#' 
#' @return The M by M item-weighted combined input matrix
#' 
#' @examples
#'data(sports)
#'np <- dim(sports)[2]
#'P <- matrix(NA,nrow=np,ncol=np)
#'P[1,] <- c(0,5,5,10,10,10,10)
#'P[2,] <- c(5,0,5,10,10,10,10)
#'P[3,] <- c(5,5,0,10,10,10,10)
#'P[4,] <- c(10,10,10,0,5,5,5)
#'P[5,] <- c(10,10,10,5,0,5,5)
#'P[6,] <- c(10,10,10,5,5,0,5)
#'P[7,] <- c(10,10,10,5,5,5,0)
#'CIW <- iwcombinpmatr(sports,w=P) 
#' 
#' @author Alessandro Albano \email{alessandro.albano@unipa.it} \cr
#' Antonella Plaia \email{antonella.plaia@unipa.it}
#' 
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28. \cr
#' Albano, A. and Plaia, A. (2021).  Element weighted Kemeny distance for  ranking data. Electronic  Journal  of  Applied Statistical Analysis, doi: 10.1285/i20705948v14n1p117
#' 
#' @seealso \code{\link{tabulaterows}} frequency distribution of a ranking data.
#' @seealso \code{\link{combinpmatr}} combined input matrix of a ranking data set.
#' 
#' @importFrom methods is
#' 
#' @export




iwcombinpmatr<-function (X,w,Wk = NULL)
{
  #if(!is.null(dim(w))){
  if(!is(dim(w),"NULL")){
    ws_k <- w
    
  } else{
    
    # original ws_k=matrix(,length(w),length(w))
    ws_k=matrix(NA,length(w),length(w))
    
    for (j in 1:length(w)) {
      for (i in 1:length(w)){
        ws_k[j,i] <- w[j]*w[i]
      }}
  }
  
  if (is.null(Wk)) {
    CI = matrix(0, ncol(X), ncol(X))
    colnames(CI) = colnames(X)
    row.names(CI) = colnames(X)
    for (i in 1:nrow(X)) {
      sm = scorematrix(t(as.matrix(X[i, ])))* ws_k
      CI = CI + sm
    }
  }
  else {
    if (is.numeric(Wk)) {
      Wk = matrix(Wk, ncol = 1)
    }
    CI = matrix(0, ncol(X), ncol(X))
    colnames(CI) = colnames(X)
    row.names(CI) = colnames(X)
    for (i in 1:nrow(X)) {
      sm = scorematrix(t(as.matrix(X[i, ])))* ws_k* Wk[i]
      CI = CI + sm
    }
  }
  CI
}
