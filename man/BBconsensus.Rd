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
BBconsensus(RR, cij)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RR}{
Candidate to be the consensus ranking
}
  \item{cij}{
Combined input matrix of the data set
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
Amodio, S., D'Ambrosio, A., & Siciliano, R. (2015). A Fast accurate heuristic algorithm for the consensus ranking problem. Submitted. Technical report available as arXiv preprint arXiv:1502.06498.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it>
}


