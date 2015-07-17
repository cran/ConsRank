\name{labels}
\alias{labels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform a ranking into a ordering. 
}
\description{
Given a ranking (or a matrix of rank data), transforms it into an ordering (or a ordering matrix)
}
\usage{
labels(x, m, label = 1:m, labs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a ranking, or a n by m data matrix in which there are n judges ranking m objects
}
  \item{m}{
the number of objects
}
  \item{label}{
optional: the name of the objects
}
  \item{labs}{
= 1 displays the names of the objects if there is argument "label", otherwise displays the permutation of first m integer.
labs = 2 is to be used only if the argument "label" is not defined. In such a case it displays the permutation of the first m letters
}
}

\value{
Ordering
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Sonia Amodio <sonia.amodio@unina.it>
}

\examples{
data(Idea)
TR=tabulaterows(Idea)
Ord=labels(TR$X,ncol(Idea),colnames(Idea),labs=1)
Ord2=labels(TR$X,ncol(Idea),labs=2)
cbind(Ord,TR$Wk)
cbind(Ord2,TR$Wk)


}

