\name{EMD}
\alias{EMD}
\docType{data}
\title{
Emond and Mason data
}
\description{
Data simuated by Emond and Mason to check their branch-and-bound algorithm. There are 112 voters ranking 15 objects. There are 21 uncomplete rankings. Data are in form of frequency table.
}
\usage{data(EMD)}


\source{
Emond, E. J., & Mason, D. W. (2000). A new technique for high level decision support. Department of National Defence, Operational Research Division, pag. 28.
}
\references{
Emond, E. J., & Mason, D. W. (2000). A new technique for high level decision support. Department of National Defence, Operational Research Division, pag. 28.
}
\examples{
data(EMD)
CR=QuickCons(EMD[,1:15],EMD[,16])

}
\keyword{datasets}
