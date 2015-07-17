\name{kemenydesign}
\alias{kemenydesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function
}
\description{
Define a design matrix to compute Kemeny distance
}
\usage{
kemenydesign(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects represented by the columns. 
}
}

\value{
Design matrix
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
D'Ambrosio, A. (2008). Tree based methods for data editing and preference rankings. Unpublished PhD Thesis. Universita' degli Studi di Napoli Federico II.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it>
}

