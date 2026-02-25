##### ITEM WEIGHTED KEMENY DISTANCE

#' Item-weighted Kemeny distance
#'
#' Compute the item-weighted Kemeny distance of a data matrix containing preference rankings, or compute the kemeny distance between two (matrices containing) rankings.
#'
#' @param x A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. If there is only x as input, the output is a square distance matrix
#' @param y A row vector, or a N by M data matrix in which there are N judges and the same M objects as x to be judged.
#' @param w A M-dimensional row vector (individually weighted items), or a M by M matrix (item similarities)
#' @return If there is only x as input, d = square distance matrix. If there is also y as input, d = matrix with N rows and n columns.
#'
#' @references Kemeny, J. G., & Snell, L. J. (1962). Preference ranking: an axiomatic approach. Mathematical models in the social sciences, 9-23. \cr
#' Albano, A. and Plaia, A. (2021)  Element weighted Kemeny distance for  ranking data. Electronic  Journal  of  Applied Statistical Analysis, doi: 10.1285/i20705948v14n1p117
#'
#' @author Alessandro Albano \email{alessandro.albano@unipa.it} \cr
#' Antonella Plaia \email{antonella.plaia@unipa.it}
#'
#' @examples
#' #Individually weighted items
#'data("German")
#'w=c(10,5,5,10)
#'iw_kemenyd(x= German[c(1,200,300,500),],w= w)
#'iw_kemenyd(x= German[1,],y=German[400,],w= w)
#'
#' #Item similarity weights
#'data(sports)
#'P=matrix(NA,nrow=7,ncol=7)
#'P[1,]=c(0,5,5,10,10,10,10)
#'P[2,]=c(5,0,5,10,10,10,10)
#'P[3,]=c(5,5,0,10,10,10,10)
#'P[4,]=c(10,10,10,0,5,5,5)
#'P[5,]=c(10,10,10,5,0,5,5)
#'P[6,]=c(10,10,10,5,5,0,5)
#'P[7,]=c(10,10,10,5,5,5,0)
#'iw_kemenyd(x=sports[c(1,3,5,7),], w= P)
#'iw_kemenyd(x=sports[1,],y=sports[100,], w= P)

#' @keywords item-weighted Kemeny distance
#'
#'
#' @importFrom proxy dist
#' @importFrom gtools combinations
#' @importFrom tidyr crossing
#' 
#' @seealso \code{\link{iw_tau_x}} item-weighted tau_x rank correlation coefficient
#' @seealso \code{\link{kemenyd}} Kemeny distance
#' 
#' @export

iw_kemenyd <- function(x,y=NULL,w){

  if(is(y,"NULL")){

    if(!is(dim(w),"NULL")){

      x <- t(apply(x,1,rank,ties.method ="min"))

      ws_k <- w
    }
    else{
      x <- x[,which(w!=0)]
      x <- t(apply(x,1,rank,ties.method ="min"))
      w <- w[which(w!=0)]
      # original ws_k <- matrix(,length(w),length(w))
      ws_k <- matrix(NA,length(w),length(w))
      for (j in 1:length(w)) {
        for (i in 1:length(w)){
          ws_k[j,i] <- w[j]*w[i]
        }}
    }

    N <- nrow(x)
    indice <- combinations(N, 2, repeats.allowed = T)
    kx <- mat.or.vec(N,N)
    for (j in 1:nrow(indice)) {
      kx[indice[j,1],indice[j,2]]= kx[indice[j,2],indice[j,1]] <- sum(abs(kemenyscore(x[ indice[j, 1],])-kemenyscore(x[ indice[j, 2],]))*ws_k/2)
    }

    kx

  }

  else {

    if (is(x, "numeric") & !is(x, "matrix")) {
      x <- matrix(x, ncol = length(x))
    }

    if (is(y, "numeric") & !is(y, "matrix")) {
      y <- matrix(y, ncol = length(y))
    }

    if(!is(dim(w),"NULL")){

      x<-t(apply(x,1,rank,ties.method ="min"))
      y<-t(apply(y,1,rank,ties.method ="min"))

      ws_k <- w }

    else{
      x <- x[,which(w!=0),drop=FALSE]
      x <- t(apply(x,1,rank,ties.method ="min"))

      y <- y[,which(w!=0),drop=FALSE]
      y <- t(apply(y,1,rank,ties.method ="min"))

      w <- w[which(w!=0)]
      # original ws_k <- matrix(,length(w),length(w))
      ws_k <- matrix(NA,length(w),length(w))

      for (j in 1:length(w)) {
        for (i in 1:length(w)){
          ws_k[j,i] <- w[j]*w[i]
        }}
    }
    nx <- nrow(x)
    ny <- nrow(y)
    indice <- data.frame(crossing(1:nx,1:ny))
    kx <- matrix(NA,nx,ny)
    for (j in 1:nrow(indice)) {
      kx[indice[j,1],indice[j,2]] <- sum(abs(kemenyscore(x[ indice[j, 1],])-kemenyscore(y[ indice[j, 2],]))*ws_k/2)
    }
    #tril(kx)
    kx
  }
}

#------------------------------------------------------------------------

# ITEM WEIGHTED TAU_X


#' Item-weighted TauX rank correlation coefficient
#'
#' Compute the item-weighted TauX rank correlation coefficient of a data matrix containing preference rankings, or compute the item-weighted correlation coefficient  between two (matrices containing) rankings.
#'
#' @param x A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. If there is only x as input, the output is a square matrix
#' @param y A row vector, or a N by M data matrix in which there are N judges and the same M objects as x to be judged.
#' @param w A M-dimensional row vector (individually weighted items), or a M by M matrix (item similarities)
#' @return Item-weighted TauX rank correlation coefficient
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28. \cr
#' Albano, A. and Plaia, A. (2021)  Element weighted Kemeny distance for  ranking data. Electronic  Journal  of  Applied Statistical Analysis, doi: 10.1285/i20705948v14n1p117
#'
#' @author Alessandro Albano \email{alessandro.albano@unipa.it} \cr
#' Antonella Plaia \email{antonella.plaia@unipa.it}
#' 
#' @seealso \code{\link{tau_x}} TauX rank correlation coefficient
#' @seealso \code{\link{iw_kemenyd}} item-weighted Kemeny distance
#'
#' @examples
#' #Individually weighted items
#'data("German")
#'w=c(10,5,5,10)
#'iw_tau_x(x= German[c(1,200,300,500),],w= w)
#'iw_tau_x(x= German[1,],y=German[400,],w= w)
#'
#' #Item similarity weights
#'data(sports)
#'P=matrix(NA,nrow=7,ncol=7)
#'P[1,]=c(0,5,5,10,10,10,10)
#'P[2,]=c(5,0,5,10,10,10,10)
#'P[3,]=c(5,5,0,10,10,10,10)
#'P[4,]=c(10,10,10,0,5,5,5)
#'P[5,]=c(10,10,10,5,0,5,5)
#'P[6,]=c(10,10,10,5,5,0,5)
#'P[7,]=c(10,10,10,5,5,5,0)
#'iw_tau_x(x=sports[c(1,3,5,7),], w= P)
#'iw_tau_x(x=sports[1,],y=sports[100,], w= P)

#' @keywords item-weighted rank correlation coefficient
#'
#' @export
#'
#' @importFrom proxy dist

iw_tau_x <- function(x,y=NULL,w){

  if(is.null(y)){

    if(!is.null(dim(w))){

      x <- t(apply(x,1,rank,ties.method ="min"))
      ws_k <- w
    }

    else{
      x <- x[,which(w!=0)]
      x <- t(apply(x,1,rank,ties.method ="min"))
      w <- w[which(w!=0)]
     #original ws_k <- matrix(,length(w),length(w))
      ws_k <- matrix(NA,length(w),length(w))

      for (j in 1:length(w)) {
        for (i in 1:length(w)){
          ws_k[j,i] <- w[j]*w[i]
        }}
    }

    N <- nrow(x)
    indice <- combinations(N, 2, repeats.allowed = T)
    tau <- mat.or.vec(N,N)
    for (j in 1:nrow(indice)) {
      tau[indice[j,1],indice[j,2]]=tau[indice[j,2],indice[j,1]] <- sum(scorematrix(x[ indice[j, 1],])*scorematrix(x[ indice[j, 2],])*ws_k)/(sum(ws_k)-sum(diag(ws_k)))

    }

    tau

  }

  else {

    if (is(x, "numeric") & !is(x, "matrix")) {
      x <- matrix(x, ncol = length(x))
    }

    if (is(y, "numeric") & !is(y, "matrix")) {
      y <- matrix(y, ncol = length(y))
    }

    if(!is.null(dim(w))){

      x <- t(apply(x,1,rank,ties.method ="min"))
      y <- t(apply(y,1,rank,ties.method ="min"))

      ws_k <- w }

    else{
      x <- x[,which(w!=0),drop=FALSE]
      x <- t(apply(x,1,rank,ties.method ="min"))

      y <- y[,which(w!=0),drop=FALSE]
      y <- t(apply(y,1,rank,ties.method ="min"))

      w <- w[which(w!=0)]
      # original ws_k <- matrix(,length(w),length(w))
      ws_k <- matrix(NA,length(w),length(w))

      for (j in 1:length(w)) {
        for (i in 1:length(w)){
          ws_k[j,i] <- w[j]*w[i]
        }}
    }
    nx <- nrow(x)
    ny <- nrow(y)
    indice <- data.frame(crossing(1:nx,1:ny))
    tau <- matrix(NA,nx,ny)
    for (j in 1:nrow(indice)) {
      tau[indice[j,1],indice[j,2]] <- sum(scorematrix(x[ indice[j, 1],])*scorematrix(y[ indice[j, 2],])*ws_k)/(sum(ws_k)-sum(diag(ws_k)))
    }

    tau
  }
}

