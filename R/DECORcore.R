#' Differential Evolution algorithm for Median Ranking
#'
#' Core function of the DECOR algorithm 
#'
#' @param cij combined input matrix
#' @param NJ the number of judjes
#' @param NP The number of population individuals
#' @param L Generations limit: maximum number of consecutive generations without improvement
#' @param FF The scaling rate for mutation. Must be in [0,1]
#' @param CR The crossover range. Must be in [0,1]
#' @param FULL Default FULL=FALSE. If FULL=TRUE, the searching is limited to the space of full rankings. In this case, the data matrix must contain full rankings.
#' 
#' @return a "list" containing the following components:
#' \tabular{lll}{
#' ConsR \tab  \tab the Consensus Ranking\cr
#' Tau \tab       \tab averaged TauX rank correlation coefficient\cr
#' besti\tab   \tab matrix of best individuals for every generation\cr
#' bestc\tab  \tab vector of best individuals' cost for every gen\cr
#' bests\tab   \tab  \cr vector of best individuals
#' avgTau\tab   \tab  \cr maximum average tauX
#' Eltime\tab  \tab Elapsed time in seconds}
#' 
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it} and Giulio Mazzeo \email{giuliomazzeo@gmail.com}
#' 
#' @references D'Ambrosio, A., Mazzeo, G., Iorio, C., and Siciliano, R. (2017). A differential evolution algorithm for finding the median ranking under the Kemeny axiomatic approach. Computers and Operations Research, vol. 82, pp. 126-138. 
#' 
#' @export




DECORcore = function(cij,NJ,NP=15,L=50,FF=0.4,CR=0.9,FULL=FALSE){
  
  # DECoR Differential Evolution for COnsensus Ranking
  
  # Discrete version of Differential Evolution algorithm
  # specifically created for solving the Consensus Ranking problem,
  # AKA social choice, AKA rank aggregation.
  # This file uses external functions for mutation and crossover during
  # the genetic evolutions of population.
  # Another function is needed to discretize and correct child solutions
  # when mutation and crossover generate out-of-bound values and duplicates.
  #
  # Input parameters:
  #
  # NP            - The number of population individuals
  #
  # L             - Generations limit: maximum number of consecutive 
  #                 generations without improvement
  #
  # F             - The scaling rate for mutation. Must be in [0,1]
  #
  # CR            - The crossover range. Must be in [0,1]
  #
  # cij           - Combined input matrix (see Emond and Mason B&B)
  #
  # NJ            - Number of judges (needed for cost computation)
  #
  # FULL          - FULL = 1 search in the space of full rankings
  #
  # Output parameters:
  #
  # besti         - Best Individuals:
  #                 matrix of best individuals for every generation
  #
  # bestc         - Best Costs:
  #                 vector of best individuals' cost for every gen
  #
  # bests         - Best Solutions:
  #                 matrix with "all" the best solutions (founded)
  #                 for checking if more than one optimal solution is found
  #
  # Notes:
  #
  # mutation      - rand/1/bin > best/1/bin
  # 
  #Authors : Giulio Mazzeo and Antonio D'Ambrosio
  
  
  # preparation
  
  tic = proc.time()[3]
  N=nrow(cij)        # number of objects                     
  costs       = matrix(0,1,NP)  # array of initial costs
  
  # initialize the population (random selection)
  population = matrix(0,(NP-1),N)
  
  for (k in 1:(NP-1)){ population[k,] = sample(N)}
  
  
  # insert a very good candidate
  
  population=rbind(population,findconsensusBB(cij))
  
  if (FULL==TRUE){
    population[NP,]=order(population[NP,])
  }
  
  
  
  # compute costs of initial population
  costs=matrix(0,NP,1)
  taos=costs
  for (i in 1:NP){
    
    COTA=combincost(population[i,],cij,NJ)
    costs[i]=COTA$cp
    taos[i]=COTA$tp
  }
  
  # store the best individual and cost of initial population
  bestc = min(costs)
  bestind = which(costs==min(costs))
  bestT = max(taos)
  besti=population[bestind,]
  
  
  # generation index
  g = 2
  
  # evolve for generations
  no_gain = 0
  
  while (no_gain < L){
    
    
    # individuals mutation
    for (i in 1:NP){
      
      
      # apply mutation
      evolution = mutaterand1(population,FF,i);
      
      # apply crossover
      evolution = crossover(population[i,],evolution,CR)
      
      
      # apply discretization and convert to rank
      
      if (FULL==TRUE){
        
        evolution=order(evolution)}
      
      else{
        
        evolution = childtie(evolution)
      }
      
      # apply selection, hold the best individual
      COTAN = combincost(evolution,cij,NJ)
      cost_new=COTAN$cp
      ta_new=COTAN$tp
      
      
      if (cost_new < costs[i]){
        population[i,] = evolution
        costs[i] = cost_new
        taos[i]=ta_new
      }
      
    }
    
    # store the best individual of current generation
    
    bestco = min(costs)
    bestc=rbind(bestc,bestco)
    bestind = which.min(costs)
    bestTa = max(taos)
    bestT=rbind(bestT,bestTa)
    bestin=population[bestind,]
    besti=rbind(besti,bestin)
    
    
    # check if this generation improved solutions
    if (bestc[g] == bestc[(g-1)]){
      
      no_gain = no_gain + 1}
    
    else{
      
      no_gain = 0
    }
    
    
    # next generation
    g = g + 1
    
  } #end while
  
  # select ALL the best solutions
  
  indexes = which(bestc==min(bestc))
  if (FULL==TRUE){ #if1
    
    if (length(indexes)==1){ #if2
      bests=childclosint(matrix(besti[indexes,],1,N))}
    else{
      bests=matrix(0,length(indexes),N)
      for (j in 1:length(indexes)){
        bests[j,]=childclosint(besti[indexes[j],])
      } #end for
    } #end if2
    
  } else { #if FULL = FALSE
    
    if(length(indexes)==1){
      
      bests = reordering(matrix(besti[indexes,],1,N))
      
    } else {
      
      bests = reordering(besti[indexes,])}
    
  } #end if1
  
  avgTau = bestT[indexes]
  
  ConsR=unique(bests)
  Tau=matrix(rep(avgTau,nrow(ConsR)),nrow(ConsR),1)
  
  
  toc = proc.time()[3]
  eltime=toc-tic
  return(list(ConsR=ConsR,Tau=Tau,besti=besti,bestc=bestc,bests=bests,avgTau=avgTau,bestT=bestT,Eltime=eltime))
  
}
