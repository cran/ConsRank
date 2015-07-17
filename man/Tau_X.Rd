\name{Tau_X}
\alias{Tau_X}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tau exstension rank correlation coefficient
}
\description{
Tau exstension is a new rank correlation coefficient defined by Emond and Mason (2002)
}
\usage{
Tau_X(X, Y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a M by N data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. 
If there is only X as input, the output is a square matrix containing the Tau_X rcc.
}
  \item{Y}{
A row vector, or a n by M data matrix in which there are n judges and the same M objects as X to be judged.
}
}

\value{
Tau_x rank correlation coefficient
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Emond, E. J., & Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.  

}
\author{
Antonio D'Ambrosio <antdambr@unina.it>
}


\seealso{
kemenyd
}
\examples{
data(BU)
RD=BU[,1:3]
Tau=Tau_X(RD)
Tau1_3=Tau_X(RD[1,],RD[3,])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Tau_x rank correlation coefficient }

