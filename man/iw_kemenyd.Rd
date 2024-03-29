% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted_KemenyDist&TauX.R
\name{iw_kemenyd}
\alias{iw_kemenyd}
\title{Item-weighted Kemeny distance}
\usage{
iw_kemenyd(x, y = NULL, w)
}
\arguments{
\item{x}{A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. If there is only x as input, the output is a square distance matrix}

\item{y}{A row vector, or a N by M data matrix in which there are N judges and the same M objects as x to be judged.}

\item{w}{A M-dimensional row vector (individually weighted items), or a M by M matrix (item similarities)}
}
\value{
If there is only x as input, d = square distance matrix. If there is also y as input, d = matrix with N rows and n columns.
}
\description{
Compute the item-weighted Kemeny distance of a data matrix containing preference rankings, or compute the kemeny distance between two (matrices containing) rankings.
}
\examples{
#Individually weighted items
data("German")
w=c(10,5,5,10)
iw_kemenyd(x= German[c(1,200,300,500),],w= w)
iw_kemenyd(x= German[1,],y=German[400,],w= w)

#Item similarity weights
data(sports)
P=matrix(NA,nrow=7,ncol=7)
P[1,]=c(0,5,5,10,10,10,10)
P[2,]=c(5,0,5,10,10,10,10)
P[3,]=c(5,5,0,10,10,10,10)
P[4,]=c(10,10,10,0,5,5,5)
P[5,]=c(10,10,10,5,0,5,5)
P[6,]=c(10,10,10,5,5,0,5)
P[7,]=c(10,10,10,5,5,5,0)
iw_kemenyd(x=sports[c(1,3,5,7),], w= P)
iw_kemenyd(x=sports[1,],y=sports[100,], w= P)
}
\references{
Kemeny, J. G., & Snell, L. J. (1962). Preference ranking: an axiomatic approach. Mathematical models in the social sciences, 9-23. \cr
Albano, A. and Plaia, A. (2021)  Element weighted Kemeny distance for  ranking data. Electronic  Journal  of  Applied Statistical Analysis, doi: 10.1285/i20705948v14n1p117
}
\seealso{
\code{\link{iw_tau_x}} item-weighted tau_x rank correlation coefficient

\code{\link{kemenyd}} Kemeny distance
}
\author{
Alessandro Albano \email{alessandro.albano@unipa.it} \cr
Antonella Plaia \email{antonella.plaia@unipa.it}
}
\keyword{Kemeny}
\keyword{distance}
\keyword{item-weighted}
