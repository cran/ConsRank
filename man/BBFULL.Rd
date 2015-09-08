\name{BBFULL}
\alias{BBFULL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Branch and Bound algorithm to find the median ranking in the space of full (or complete) rankings.
}
\description{
Branch-and-bound algorithm to find consensus ranking as definned by D'Ambrosio et al. (2015). If the number of objects to be ranked is large (greater than 20 or 25), it can work for very long time. Use either QuickCons or FASTcons with the option FULL=TRUE instead
}
\usage{
BBFULL(X,Wk=NULL,PS=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A N by M data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. The data matrix MUST contain only full rankings, or incomplete rankings, without ties.\\
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
\item{Tau }{Averaged Kendall's Tau}
\item{Eltime }{Elapsed time in seconds}
%% ...
}
\references{
D'Ambrosio, A., Amodio, S. & Iorio, C. (2015). Two algorithms for finding optimal solutions of the Kemeny rank
aggregation problem for full rankings. Technical report, University of Naples Federico II.
}
\author{
Antonio D'Ambrosio <antdambr@unina.it> 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
FASTcons

QuickCons
}
\examples{
data(APAFULL)
CR=BBFULL(APAFULL)

}
\keyword{ Consensus ranking }
\keyword{ Median ranking }% __ONLY ONE__ keyword per line
