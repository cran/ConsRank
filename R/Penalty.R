#' Auxiliary function
#'
#' Assign a penalty to the branches of the FASTcons amd QuickCons algorithms 
#'
#' @param CR candidate to be the median ranking
#' @param cij combined input matrix
#' @param indice other input called by other functions
#'
#' @return a penalty value
#'
#' @references Amodio, S., D'Ambrosio, A. and Siciliano, R. (2016). Accurate algorithms for identifying the median ranking when dealing with weak and partial rankings under the Kemeny axiomatic approach. European Journal of Operational Research, 249(2), 667-676.
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Sonia Amodio \email{sonia.amodio@unina.it}
#' 
#' @export



Penalty = function(CR,cij,indice)   #indice must be order(CR)
{
  if (CR[indice[1,1]] < CR[indice[1,2]]) { #case 1, the first object is preferred
    #to the second object
    if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
      Po=0
    } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po = cij[indice[1,2],indice[1,1]]-cij[indice[1,1],indice[1,2]]  #     cji-cij
    } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po = cij[indice[1,2],indice[1,1]] #cji
    }
  } else if (CR[indice[1,1]] > CR[indice[1,2]]) { #case 2 the first object is not
    #preferred to the second one
    if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
      Po= cij[indice[1,1],indice[1,2]]-cij[indice[1,2],indice[1,1]]  #cij-cji
    } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po = 0
    } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po = cij[indice[1,1],indice[1,2]] #cij
    }
  } else if (CR[indice[1,1]] == CR[indice[1,2]]) { #case 3 they are in a tie
    if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == -1) {
      Po = -cij[indice[1,2],indice[1,1]] #-cj
    } else if (sign(cij[indice[1,1],indice[1,2]]) == -1 & sign (cij[indice[1,2],indice[1,1]]) == 1) {
      Po = -cij[indice[1,1],indice[1,2]] #-cij
    } else if (sign(cij[indice[1,1],indice[1,2]]) == 1 & sign (cij[indice[1,2],indice[1,1]]) == 1 |
               sign(cij[indice[1,1],indice[1,2]]) == 0 & sign (cij[indice[1,2],indice[1,1]]) == 0 )  {
      Po = 0
    }
    
  }
  Po
}

