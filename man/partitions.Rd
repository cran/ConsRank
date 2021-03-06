% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partitions.R
\name{partitions}
\alias{partitions}
\title{Generate partitions of n items constrained into k non empty subsets}
\usage{
partitions(n, k = NULL, items = NULL, itemtype = "L")
}
\arguments{
\item{n}{a (integer) number denoting the number of items}

\item{k}{The number of the non-empty subsets. Default value is NULL, in this case all the possible partitions are displayed}

\item{items}{items: the items to be placed into the ordering matrix. Default are the first c small letters}

\item{itemtype}{to be used only if items is not set. The default value is "L", namely letters. Any other symbol produces items as the first c integers}
}
\value{
the ordering matrix (or vector)
}
\description{
Generate all possible partitions of n items constrained into k non empty subsets. It does not generate the universe of rankings constrained into k buckets.
}
\details{
If the objects to be ranked is large (>15-20) with some missing, it can take long time to find the solutions. If the searching space is 
limited to the space of full rankings (also incomplete rankings, but without ties), use the function BBFULL or the functions FASTcons and QuickCons 
with the option FULL=TRUE.
}
\examples{
X<-partitions(4,3)
#shows all the ways to partition 4 items (say "a", "b", "c" and "d" into 3 non-empty subets
 #(i.e., into 3 buckets). The Stirling number of the second kind (4,3) indicates that there
 #are 6 ways.
s2<-stirling2(4,3)$S
X2<-order2rank(X) #it transform the ordering into ranking


}
\seealso{
\code{\link{stirling2}} Stirling number of second kind.

\code{\link{rank2order}} Convert rankings into orderings.

\code{\link{order2rank}} Convert orderings into ranks.

\code{\link{univranks}} Generate the universe of rankings given the input partition
}
\author{
Antonio D'Ambrosio \email{antdambr@unina.it}
}
