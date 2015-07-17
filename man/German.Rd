\name{German}
\alias{German}
\docType{data}
\title{
German political goals
}
\description{
Ranking data of 2262 German respondents about the desirability of the four political goals: a = the maintenance of order in the nation; b = giving people more say in the decisions of government; c = growthing rising prices; d = protecting freedom of speech
}
\usage{data(German)}


\source{
Croon, M. A. (1989). Latent class models for the analysis of rankings. Advances in psychology, 60, 99-121.
}

\examples{
data(German)
TR=tabulaterows(German)
polyplot(TR$X,Wk=TR$Wk,nobj=4)
}
\keyword{datasets}
