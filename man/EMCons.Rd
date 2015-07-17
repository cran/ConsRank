\name{EMCons}
\alias{EMCons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Branch-and-bound algorithm to find consensus (median) ranking according to the Kemeny's axiomatic approach
}
\description{
Branch-and-bound algorithm to find consensus ranking as definned by Emond and Mason (2002). If the number of objects to be ranked is large (greater than 15 or 20), it can work for very long time.
}
\usage{
EMCons(X, Wk = NULL, PS = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. \\
Alternatively X can contain the rankings observed only once. In this case the argument Wk must be used
}
  \item{Wk}{
Optional: the frequency of each ranking in the data
}
  \item{PS}{
If PS=TRUE, on the screen some information about how many branches are processed are displayed
}
}
\details{
If the objects to be ranked is large (>15-20), it can take long time to finde the solutions
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Consensus }{Consensus ranking}
\item{Tau }{Tau extended rank correlation coefficient}
\item{Eltime }{Elapsed time in seconds}
%% ...
}
\references{
Emond, E. J., & Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it> and Sonia Amodio <sonia.amodio@unina.it>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
FASTcons

QuickCons
}
\examples{
data(Idea)
RevIdea=6-Idea
# as 5 means "most associated", it is necessary compute the reverse ranking 
#of each rankings to have rank 1 = "most associated" and rank 5 = "least associated"
CR=EMCons(RevIdea)

}
\keyword{ Consensus ranking }
\keyword{ Median ranking }% __ONLY ONE__ keyword per line
