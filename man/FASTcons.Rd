\name{FASTcons}
\alias{FASTcons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FAST algorithm to find consensus (median) ranking according the Kemeny's axionatic approach
}
\description{
FAST algorithm to find consensus (median) ranking defined by Amodio, D'Ambrosio and Siciliano (2015). It returns at least one of the solutions. If there are multiple solutions, sometimes it returns all the solutions, sometimes it returns some solutions, always it returns at least one solution.
}
\usage{
FASTcons(X, Wk = NULL, maxiter = 50, FULL=FALSE)
}

\arguments{
  \item{X}{
A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. 
Alternatively X can contain the rankings observed only once in the sample. In this case the argument Wk must be used
}
  \item{Wk}{
Optional: the frequency of each ranking in the data
}
  \item{maxiter}{
maximum number of iterations: default = 50.
}
  \item{FULL}{
Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings. In this case, the data matrix must contain full rankings. 
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Consensus }{Consensus ranking}
\item{Tau }{Tau extended rank correlation coefficient}
\item{Eltime }{Elapsed time in seconds}
%% ...
}
\references{
Amodio, S., D'Ambrosio, A. & Siciliano, R (2015). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research. DOI: 10.1016/j.ejor.2015.08.048.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it> and Sonia Amodio <sonia.amodio@unina.it>
}



\seealso{
EMCons

QuickCons


}
\examples{
#data(EMD)
#X=EMD[,1:15]
#Wk=matrix(EMD[,16],nrow=nrow(X))
#CR=FASTcons(X,Wk,maxiter=100)
#These lines produce all the three solutions in less than a minute.

data(sports)
CR=FASTcons(sports,maxiter=10)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FAST algorithm }

