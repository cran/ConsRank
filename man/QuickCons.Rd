\name{QuickCons}
\alias{QuickCons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quick algorithm to find up to 4 solutions to the consensus ranking problem
}
\description{
The Quick algorithm finds up to 4 solutions. Solutions reached are most of the time optimal solutions. 
}
\usage{
QuickCons(X, Wk = NULL, FULL = FALSE, PS=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A N by M data matrix in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns.
Alternatively X can contain the rankings observed only once in the sample. In this case the argument Wk must be used
}
  \item{Wk}{
Optional: the frequency of each ranking in the data
}
  \item{FULL}{
Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings. In this case, the data matrix must contain full rankings. 
}
  \item{PS}{
Default PS=FALSE. If PS=TRUE the number of evaluated branches is diplayed
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Consensus }{Consensus ranking}
\item{Tau }{Tau extension rank correlation coefficient}
%% ...
}
\references{
Amodio, S., D'Ambrosio, A. & Siciliano, R (2015). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research. DOI: 10.1016/j.ejor.2015.08.048.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it> and Sonia Amodio <sonia.amodio@unina.it>
}




\seealso{
FASTcons

EMCons
}
\examples{
data(EMD)
CR=QuickCons(EMD[,1:15],EMD[,16])



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Quick algorithm }

