\name{ConsRank-package}
\alias{ConsRank-package}
\alias{ConsRank}
\docType{package}
\title{
Compute the Median Ranking(s) According to the Kemeny's Axiomatic Approach
}
\description{
Compute the median ranking according the Kemeny's axiomatic approach. Rankings can or cannot contain ties, rankings can be both complete or incomplete. 
}
\details{
\tabular{ll}{
Package: \tab ConsRank\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2015-09-01\cr
License: \tab GPL-3\cr
}

}
\author{
Antonio D'Ambrosio <antdambr@unina.it>, Sonia Amdio <sonia.amodio@unina.it>

Maintainer: Antonio D'Ambrosio <antdambr@unina.it>
}
\references{
Kemeny, J. G., & Snell, J. L. (1962). Mathematical models in the social sciences (Vol. 9). New York: Ginn.

Marden, J. I. (1996). Analyzing and modeling rank data. CRC Press.

Emond, E. J., & Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.

Heiser, W. J., & D'Ambrosio, A. (2013). Clustering and prediction of rankings within a Kemeny distance framework. In Algorithms from and for Nature and Life (pp. 19-31). Springer International Publishing.

Amodio, S., D'Ambrosio, A. & Siciliano, R (2015). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research. DOI: 10.1016/j.ejor.2015.08.048.

D'Ambrosio, A., Amodio, S. & Iorio, C. (2015). Two algorithms for finding optimal solutions of the Kemeny rank
aggregation problem for full rankings. Technical report, University of Naples Federico II. Submitted.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%~~ documentation directory ~~
\keyword{ package }

\examples{
## load APA data set, full version
data(APAFULL)
## Emond and Mason Branch-and-Bound algorithm. 
## If the number of object is higher than 20, EMCons function may work for several minutes. 
## Use either QuickCons or FASTcons instead
CR=EMCons(APAFULL)
TR=tabulaterows(APAFULL)
CR2=FASTcons(TR$X,TR$Wk,maxiter=5)
CR3=QuickCons(TR$X,TR$Wk)

#####################################
### load sports data set
#data(sports)
### FAST algorithm
#CR=FASTcons(sports,maxiter=10)
#####################################

#######################################
### load Emond and Mason data set
#data(EMD)
### matrix X contains rankings
#X=EMD[,1:15]
### vector Wk contains frequencies
#Wk=EMD[,16]
### QUICK algorithm
#CR=QuickCons(X,Wk)
#######################################

}
