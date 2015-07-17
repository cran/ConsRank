\name{BU}
\alias{BU}
\docType{data}
\title{
Brook and Upton data
}
\description{
The data consist of ballots of three candidates, where the 948 voters rank the candidates from 1 to 3. Data are in form of frequency table.
}
\usage{data(BU)}

\source{
Brook, D., & Upton, G. J. G. (1974). Biases in local government elections due to position on the ballot paper. Applied Statistics, 414-419.
}
\references{
Marden, J. I. (1996). Analyzing and modeling rank data. CRC Press, pag. 153.
}
\examples{
data(BU)
polyplot(BU[,1:3],Wk=BU[,4])

}
\keyword{datasets}
