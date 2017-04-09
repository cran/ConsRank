#' Frequency distribution of a sample of rankings
#'
#' Given a sample of preference rankings, it compute the frequency associated to each ranking
#'
#' @param X a N by M data matrix containing N judges judging M objects
#' @param miss TRUE if there are missing data (either partial or incomplete rankings): default: FALSE
#'
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' X \tab  \tab the unique rankings\cr
#' Wk \tab       \tab the frequency associated to each ranking\cr
#' tabfreq\tab   \tab frequency table}
#' 
#' @examples
#' data(Idea)
#' TR=tabulaterows(Idea)
#' FR=TR$Wk/sum(TR$Wk)
#' RF=cbind(TR$X,FR)
#' colnames(RF)=c(colnames(Idea),"fi")
#' #compute modal ranking
#' maxfreq=which(RF[,6]==max(RF[,6]))
#' labels(RF[maxfreq,1:5],5,colnames(Idea),labs=1)
#' #
#' data(APAred)
#' TR=tabulaterows(APAred)
#' #
#' data(APAFULL)
#' TR=tabulaterows(APAFULL)
#' CR1=EMCons(TR$X,TR$Wk)
#' CR2=FASTcons(TR$X,TR$Wk,maxiter=15)
#' CR3=QuickCons(TR$X,TR$Wk)
#' 
#' @author Sonia Amodio \email{sonia.amodio@unina.it}
#' 
#' @keywords modal ranking
#' 
#' @export




tabulaterows = function(X,miss=FALSE) {
  
  #given a sample of preference rankings, it counts the judges that have equal preferences
  #and it tabulates the row of the data matrix
  
  if (sum(is.na(X))>0) {
    miss=TRUE
    X[is.na(X)]=-10
  }
  
  coun = table(apply(X, 1, paste, collapse=","))
  nam = names(coun)
  spl = (strsplit(nam, ","))
  kkn  = lapply(spl, as.numeric)
  tab = t(as.data.frame(kkn))
  cek =  cbind(tab,coun)
  coun=as.matrix(coun)
  rownames(coun)=NULL
  rownames(tab)=NULL
  if (miss==TRUE) {
    tab[tab==-10]=NA
  }
  
  
  return(list(X=tab, Wk=coun, tabfreq=cbind(tab,coun)))
}

