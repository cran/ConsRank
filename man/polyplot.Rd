\name{polyplot}
\alias{polyplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot rankings on a permutation polytope of 3 o 4 objects containing all possible ties
}
\description{
Plot rankings a permutation polytope that is the geometrical space of preference rankings. The plot is available for 3 or for 4 objects
}
\usage{
polyplot(X = NULL, L = NULL, Wk = NULL, nobj = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the sample of rankings, most of the time as returned by tabulaterows
}
  \item{L}{
labels of the objects
}
  \item{Wk}{
frequency associated to each ranking
}
  \item{nobj}{
number of objects: 3 or 4
}
}
\details{
polyplot() plots the universe of 3 objecys. polyplot(nobj=4) plots the universe of 4 objecys. 
}

\references{
Thompson, G. L. (1993). Generalized permutation polytopes and exploratory graphical methods for ranked data. The Annals of Statistics, 1401-1430.

Heiser, W. J., & D'Ambrosio, A. (2013). Clustering and prediction of rankings within a Kemeny distance framework. In Algorithms from and for Nature and Life (pp. 19-31). Springer International Publishing.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it> and Sonia Amodio <sonia.amodio@unina.it>
}


\seealso{
tabulaterows
}
\examples{
polyplot()

#polyplot(nobj=4)

data(BU)

polyplot(BU[,1:3],Wk=BU[,4])


}
\keyword{ Permutation polytopes }

