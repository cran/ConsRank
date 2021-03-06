% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scorematrix.R
\name{scorematrix}
\alias{scorematrix}
\title{Score matrix according Emond and Mason (2002)}
\usage{
scorematrix(X)
}
\arguments{
\item{X}{a ranking (must be a row vector or, better, a matrix with one row and M columns)}
}
\value{
the M by M score matrix
}
\description{
Given a ranking, it computes the score matrix as defined by Emond and Mason (2002)
}
\examples{
Y <- matrix(c(1,3,5,4,2),1,5)
SM<-scorematrix(Y)
#
Z<-c(1,2,4,3)
SM2<-scorematrix(Z)

}
\references{
Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
}
\seealso{
\code{\link{combinpmatr}} The combined inut matrix
}
\author{
Antonio D'Ambrosio \email{antdambr@unina.it}
}
