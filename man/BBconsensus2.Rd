\name{BBconsensus2}
\alias{BBconsensus2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Core function in computing consensus ranking as defined by Emond and Mason (2002)
}
\description{
Core function in computing consensus ranking as defined by Emond and Mason (2002), recalled by EMCons function
}
\usage{
BBconsensus2(RR, cij, Po, PS = TRUE,FULL=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RR}{
A ranking
}
  \item{cij}{
Combined input matrix
}
  \item{Po}{
Current penalty
}
  \item{PS}{
If true, it prints the evaluating branches
}
  \item{FULL}{
Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. FULL=TRUE if the function is called by BBFULL algorithm. 
}
}

\value{
Consensus ranking
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Emond, E. J., & Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it> and Sonia Amodio <sonia.amodio@unina.it>
}

\seealso{
EMCons
}
