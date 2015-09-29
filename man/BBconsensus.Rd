\name{BBconsensus}
\alias{BBconsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the first approximation to the consensus ranking. Most of the time the output
is a solution, maybe not unique
}
\description{
Find a first approximation to the consensus ranking.
}
\usage{
BBconsensus(RR, cij, FULL=FALSE, PS=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RR}{
Candidate to be the consensus ranking
}
  \item{cij}{
Combined input matrix of the data set
}
  \item{FULL}{
Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. FULL=TRUE if the function is called by BBFULL algorithm. 
}
  \item{PS}{
Default PS=FALSE. If PS=TRUE the number of evaluated branches is diplayed
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{cons }{a first approximation to be the consensus ranking}
\item{pen }{penalty value}
%% ...
}
\references{
Amodio, S., D'Ambrosio, A. & Siciliano, R (2015). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research. DOI: 10.1016/j.ejor.2015.08.048.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it>
}


