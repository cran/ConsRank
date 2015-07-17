\name{Idea}
\alias{Idea}
\docType{data}
\title{
Idea data set
}
\description{
 98 college students where asked to rank five words, (thought,  play,  theory, dream,  attention) regarding its association with the word idea, from 5=most associated to 1=least associated.
}
\usage{data(Idea)}

\source{
Fligner, M. A., & Verducci, J. S. (1986). Distance based ranking models. Journal of the Royal Statistical Society. Series B (Methodological), 359-369.
}

\examples{
data(Idea)
revIdea=6-Idea
TR=tabulaterows(revIdea)
CR=QuickCons(TR$X,TR$Wk)
colnames(CR$Consensus)=colnames(Idea)

}
\keyword{datasets}
