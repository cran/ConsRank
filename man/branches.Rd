\name{branches}
\alias{branches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxliary code recalled by other routines
}
\description{
Branches discovery 
}
\usage{
branches(brR, cij, b, Po, ord, Pb, FULL=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{brR}{
Current processed branche of the BB algorithm
}
  \item{cij}{
Combined input matrix
}
  \item{b}{
Other inputs recalled by main functions
}
  \item{Po}{
Other inputs recalled by main functions
}

  \item{ord}{
Other inputs recalled by main functions
}

  \item{Pb}{
Other inputs recalled by main functions
}
  \item{FULL}{
Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of complete rankings. TRUE=TRUE if the function is called by BBFULL algorithm. 
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{cr}{ranking belongs to the branche}
\item{pcr }{penalty of the current ranking}
%% ...
}

\author{
Antonio D'Ambrosio <antdambr@unina.it>
}