\name{combinpmatr}
\alias{combinpmatr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combined input matrix of a data set
}
\description{
Compute the Combined input matrix of a data set as defined by Emond and Mason (2002)
}
\usage{
combinpmatr(X, Wk = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix N by M, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. \\
Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
}
  \item{Wk}{
Optional: the frequency of each ranking in the data
}
}

\value{
Combined input matrix
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
tabulaterows
}
\examples{
data(APAred) 
CI=combinpmatr(APAred) 
TR=tabulaterows(APAred) 
CI=combinpmatr(TR$X,TR$Wk)

}