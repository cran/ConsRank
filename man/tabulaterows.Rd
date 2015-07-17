\name{tabulaterows}
\alias{tabulaterows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Frequency distribution of a sample of rankings
}
\description{
Given a sample of preference rankings, it compute the frequency associated to each ranking
}
\usage{
tabulaterows(X, miss=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a N by M data matrix containing N judges judging M objects}
\item{miss}{TRUE if there are missing data (partial or uncomplete rankings): default: FALSE}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{X }{the unique rankings }
\item{Wk }{the frequency associated to each ranking}
\item{tabfreq }{frequency table}
%% ...
}

\author{
Sonia Amodio <sonia.amodio@unina.it>
}

\examples{
data(Idea)
TR=tabulaterows(Idea)
FR=TR$Wk/sum(TR$Wk)
RF=cbind(TR$X,FR)
colnames(RF)=c(colnames(Idea),"fi")
#compute modal ranking
maxfreq=which(RF[,6]==max(RF[,6]))
labels(RF[maxfreq,1:5],5,colnames(Idea),labs=1)


data(APAred)
TR=tabulaterows(APAred)

data(APAFULL)
TR=tabulaterows(APAFULL)
CR1=EMCons(TR$X,TR$Wk)
CR2=FASTcons(TR$X,TR$Wk,maxiter=15)
CR3=QuickCons(TR$X,TR$Wk)


}
\keyword{ modal ranking }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
