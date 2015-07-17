\name{kemenyd}
\alias{kemenyd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kemeny distance
}
\description{
Compute the Kemeny distance of a data matrix containing preference rankings, or compute the kemeny distance between two (matrices containing) rankings. 
}
\usage{
kemenyd(X, Y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. 
If there is only X as input, the output is a square distance matrix

}
  \item{Y}{
A row vector, or a n by M data matrix in which there are n judges and the same M objects as X to be judged.
}
}

\value{
If there is anly X as input, d = square distance matrix.

If there is also Y as input, d = matrix with N rows and n columns.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kemeny, J. G., & Snell, L. J. (1962). Preference ranking: an axiomatic approach. Mathematical models in the social sciences, 9-23.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Tau_X
}
\examples{
data(Idea)
KD=kemenyd(Idea)

KD2=kemenyd(Idea[1:10,],Idea[55,])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kemeny distance }

