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
QuickCons(X, Wk = NULL)
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
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Consensus }{Consensus ranking}
\item{Tau }{Tau extension rank correlation coefficient}
%% ...
}
\references{
Amodio, S., D'Ambrosio, A., & Siciliano, R. (2015). A Fast accurate heuristic algorithm for the consensus ranking problem. Submitted. Technical report available as arXiv preprint arXiv:1502.06498.
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

