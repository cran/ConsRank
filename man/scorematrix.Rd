\name{scorematrix}
\alias{scorematrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Score matrix according Emond and Mason (2002)
}
\description{
Given a ranking, it computes the score matrix as defined by Emond and Mason (2002)
}
\usage{
scorematrix(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a ranking (must be a row vector or, better, a matrix with one row and M columns)
}
}

\value{
the M by M score matrix
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

\examples{
Y = matrix(c(1,3,5,4,2),1,5)
SM=scorematrix(Y)

Z=c(1,2,4,3)
SM2=scorematrix(Z)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Score matrix }

