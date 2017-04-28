\name{ConsRank-package}
\alias{ConsRank-package}
\alias{ConsRank}
\docType{package}
\title{
Compute the Median Ranking According to the Kemeny's Axiomatic Approach
}
\description{
Compute the median ranking according to the Kemeny's axiomatic approach. Rankings can or cannot contain ties, rankings can be both complete or incomplete. The package contains both branch-and-bound and heuristic solutions
}
\details{
\tabular{ll}{
Package: \tab ConsRank\cr
Type: \tab Package\cr
Version: \tab 2.0.1\cr
Date: \tab 2017-04-28\cr
License: \tab GPL-3\cr
}

}
\author{
Antonio D'Ambrosio <antdambr@unina.it>, Sonia Amdio <sonia.amodio@unina.it>, Giulio Mazzeo <giuliomazzeo@gmail.com>

Maintainer: Antonio D'Ambrosio <antdambr@unina.it>
}
\references{
Kemeny, J. G., & Snell, J. L. (1962). Mathematical models in the social sciences (Vol. 9). New York: Ginn.

Marden, J. I. (1996). Analyzing and modeling rank data. CRC Press.

Emond, E. J., & Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.

D'Ambrosio, A. (2008). Tree based methods for data editing and preference rankings. Ph.D. thesis. \url{http://www.fedoa.unina.it/id/eprint/2746}

Heiser, W. J., & D'Ambrosio, A. (2013). Clustering and prediction of rankings within a Kemeny distance framework. In Algorithms from and for Nature and Life (pp. 19-31). Springer International Publishing.

Amodio, S., D'Ambrosio, A. & Siciliano, R (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, vol. 249(2).

D'Ambrosio, A., Amodio, S. & Iorio, C. (2015). Two algorithms for finding optimal solutions of the Kemeny rank aggregation problem for full rankings. Electronic Journal of Applied Statistical Analysis, vol. 8(2).

D'Ambrosio, A., Mazzeo, G., Iorio, C., & Siciliano, R. (2017). A differential evolution algorithm for finding the median ranking under the Kemeny axiomatic approach. Computers & Operations Research, vol. 82.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%~~ documentation directory ~~
\keyword{ Median ranking }
\keyword{ Consensus ranking }
\keyword{ Kemeny distance }
\keyword{ Tau_X rank correlation coefficient }
\keyword{ Branch-and-bound algorithms }
\keyword{ Differential evolution algorithms }

\examples{
## load APA data set, full version
data(APAFULL)
## Emond and Mason Branch-and-Bound algorithm. 
## If the number of object is higher than 20, EMCons function may work for several minutes. 
## Use either QuickCons, DECOR, FASTcons or FASTDECOR instead
CR=EMCons(APAFULL)
TR=tabulaterows(APAFULL)
CR2=FASTcons(TR$X,TR$Wk,maxiter=5)
CR3=QuickCons(TR$X,TR$Wk)
CR4=DECOR(TR$X,TR$Wk)

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
