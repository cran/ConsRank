\name{findbranches}
\alias{findbranches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function
}
\description{
Find correct branches in the Branch-and-Bound algorithms
}
\usage{
findbranches(R, ord, b, FULL=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
Candidate to be the consensus ranking
}
  \item{ord}{
 other input values recalled by other routines

}

  \item{b}{
 other input values recalled by other routines

}
  \item{FULL}{
Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. TRUE=TRUE if the function is called by BBFULL algorithm. 
}
}
\value{
Candidate to be the consensus ranking
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Antonio D'Ambrosio <antdambr@unina.it> and Sonia amodio sonia.amodio@unina.it>
}


