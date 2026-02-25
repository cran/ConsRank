#' Median Constrained Bucket Order (MCBO)
#'
#' Find the median ranking constrained to exactly b buckets (tied groups) according 
#' to Kemeny's axiomatic approach. This implements the algorithms described in 
#' D'Ambrosio et al. (2019) for rank aggregation with fixed bucket constraints.
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. 
#'   Each row is a ranking of the objects which are represented by the columns. 
#'   If X contains the rankings observed only once, the argument wk can be used
#' @param nbuckets Integer. The number of buckets (tied groups) the consensus ranking 
#'   must contain. Must be between 2 and M-1 (where M is the number of objects)
#' @param wk Optional: the frequency of each ranking in the data
#' @param ps Logical. If TRUE, displays progress information on screen. Default TRUE
#' @param algorithm Character string specifying the algorithm to use. One of:
#'   \itemize{
#'     \item "BB" - Branch-and-Bound (exact, default). Best for M <= 15
#'     \item "quick" - Quick algorithm (heuristic). Best for 15 < M <= 50
#'     \item "decor" - Differential Evolution (metaheuristic). Best for M > 50
#'   }
#' @param itermax Integer. Maximum number of iterations for "quick" and "decor" algorithms. 
#'   Default 10
#' @param np Integer. Number of population individuals for "decor" algorithm. Default 10
#' @param gl Integer. Generations limit for "decor": maximum number of consecutive 
#'   generations without improvement. Default 100
#' @param ff Numeric. Scaling rate for mutation in "decor". Must be in [0,1]. Default 0.4
#' @param cr Numeric. Crossover range for "decor". Must be in [0,1]. Default 0.8
#' @param use_cpp Logical. If TRUE (default), use optimized C++ implementations 
#'   for core functions (combinpmatr, scorematrix, PenaltyBB2, etc.)
#'
#' @return A list containing the following components:
#' \tabular{lll}{
#'   Consensus \tab \tab The consensus ranking with exactly nbuckets buckets\cr
#'   Tau \tab \tab Averaged TauX rank correlation coefficient\cr
#'   Eltime \tab \tab Elapsed time in seconds
#' }
#'
#' @details 
#' The median constrained bucket order problem finds the ranking that minimizes the 
#' sum of Kemeny distances to the input rankings, subject to the constraint that 
#' the solution must have exactly \code{nbuckets} groups of tied items.
#' 
#' This is useful in applications where the output must conform to predetermined 
#' categories, such as:
#' \itemize{
#'   \item Wine quality classifications (e.g., 5 fixed tiers)
#'   \item Medical triage systems (fixed severity codes)
#'   \item Educational grading (fixed letter grades: A, B, C, D, F)
#' }
#' 
#' The search space is restricted to \eqn{Z^{n/b}}, which contains
#' \deqn{\sum_{b=1}^{n} b! S(n,b)}
#' rankings, where S(n,b) is the Stirling number of the second kind.
#' 
#' \strong{Algorithm Selection Guidelines:}
#' \itemize{
#'   \item \strong{BB}: Exact solution, guaranteed optimal. Use for M <= 15 items
#'   \item \strong{quick}: Fast heuristic, near-optimal. Use for 15 < M <= 50 items
#'   \item \strong{decor}: Best for large problems. Use for M > 50 items
#' }
#' 
#' For stochastic algorithms (quick, decor), consider running multiple times 
#' (controlled by \code{itermax}) to avoid local optima.
#'
#' @examples
#' \dontrun{
#' # Simple example with 98 judges ranking 5 items into 3 buckets
#' data(Idea)
#' RevIdea <- 6 - Idea  # Reverse ranking
#' CR <- mcbo(RevIdea, nbuckets = 3, algorithm = "BB")
#' print(CR$Consensus)
#' print(CR$Tau)
#' 
#' # Large dataset with Quick algorithm
#' data(EMD)
#' CR_quick <- mcbo(EMD[,1:15], nbuckets = 5, wk = EMD[,16], 
#'                  algorithm = "quick", itermax = 20)
#' }
#'
#' @references 
#' D'Ambrosio, A., Iorio, C., Staiano, M., and Siciliano, R. (2019). 
#' Median constrained bucket order rank aggregation. 
#' Computational Statistics, 34(2), 787-802. 
#' \doi{10.1007/s00180-018-0858-z}
#' 
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso \code{\link{consrank}} for unconstrained consensus ranking
#' @seealso \code{\link{combinpmatr}} for computing the combined input matrix
#' @seealso \code{\link{scorematrix}} for computing score matrices
#' @seealso \code{\link{tau_x}} for TauX correlation coefficient
#' @seealso \code{\link{kemenyd}} for Kemeny distance
#' @seealso \code{\link{stirling2}} for Stirling numbers (bucket combinatorics)
#'
#' @keywords Median constrained bucket order
#' @keywords Median ranking
#' @keywords Consensus ranking
#' @keywords Bucket order
#' @keywords Kemeny distance
#' @keywords Rank aggregation
#'
#' @export

mcbo <- function(X, nbuckets, wk = NULL, ps = TRUE, 
                 algorithm = "BB", itermax = 10,
                 np = 10, gl = 100, ff = 0.4, cr = 0.8,
                 use_cpp = TRUE) {
  
  # ══════════════════════════════════════════════════════════════════════════
  # 1. INPUT VALIDATION
  # ══════════════════════════════════════════════════════════════════════════
  
  if (is(X, "data.frame")) X <- as.matrix(X)
  if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
  
  M <- nrow(X)
  N <- ncol(X)
  
  # Validation checks
  if (M < 2) stop("Data matrix must contain at least two different rankings")
  if (nbuckets < 2) stop("nbuckets must be at least 2")
  if (nbuckets >= N) stop(paste0("nbuckets must be less than N (", N, ")"))
  
  valid_algorithms <- c("BB", "quick", "decor")
  if (!algorithm %in% valid_algorithms) {
    stop(paste0("algorithm must be one of: ", paste(valid_algorithms, collapse = ", ")))
  }
  
  # Handle weights
  if (!is.null(wk)) {
    if (is.numeric(wk)) wk <- matrix(wk, ncol = 1)
    if (length(wk) != M) stop("Length of wk must equal number of rows in X")
    if (any(wk <= 0)) stop("All weights must be positive")
  }
  
  # ══════════════════════════════════════════════════════════════════════════
  # 2. COMPUTE COMBINED INPUT MATRIX
  # ══════════════════════════════════════════════════════════════════════════
  
  tic <- proc.time()[3]
  
  cij <- combinpmatr(X, Wk = wk, use_cpp = use_cpp)
  nj <- if (is.null(wk)) M else sum(wk)
  
  # Trivial case: all zeros
  if (sum(cij == 0) == length(cij)) {
    cr <- matrix(rep(1:nbuckets, length.out = N), nrow = 1)
    colnames(cr) <- colnames(X)
    return(list(Consensus = cr, Tau = 0, Eltime = 0))
  }
  
  # ══════════════════════════════════════════════════════════════════════════
  # 3. ALGORITHM DISPATCH
  # ══════════════════════════════════════════════════════════════════════════
  
  if (algorithm == "BB") {
    out <- consrankBBbuckets(cij, nbuckets, nj, wk, ps, use_cpp)
    
  } else if (algorithm == "quick") {
    out <- Quick_buckets(cij = cij, 
                         nbuckets = nbuckets, 
                         nj = nj, 
                         Wk = wk, 
                         maxiter = itermax,
                         PS = ps,
                         use_cpp = use_cpp)
    
  } else if (algorithm == "decor") {
    out <- DECoR_buckets(
      cij = cij,
      nbuckets = nbuckets,
      nj = nj,
      Wk = wk,
      np = np,
      gl = gl,
      ff = ff,
      cr = cr,
      maxiter = itermax,
      PS = ps,
      use_cpp = use_cpp
    )
  }
  
  toc <- proc.time()[3]
  out$Eltime <- toc - tic
  
  return(out)
}


# ══════════════════════════════════════════════════════════════════════════
# BB ALGORITHM - ORIGINAL LOGIC
# ══════════════════════════════════════════════════════════════════════════

consrankBBbuckets <- function(cij, buckets, nj, Wk, ps, use_cpp) {
  
  N <- ncol(cij)
  
  # Generate initial candidates
  R_init_res <- N2R(N, buckets, cij, nj, repts = 50, use_cpp = use_cpp)
  R_init <- matrix(R_init_res$R[1, ], nrow = 1)
  
  # Phase 1: Find initial approximation
  CA <- BBconsensusBuckts(R_init, cij, ps, buckets, use_cpp)
  Po <- CA$pen
  
  # Phase 2: Refine with full BB
  consensus <- BBconsensus3(R_init, cij, Po, ps, buckets, use_cpp)
  
  # Calculate Tau
  if (nrow(consensus) > 1) {
    TauX <- numeric(nrow(consensus))
    for (k in 1:nrow(consensus)) {
      Sij <- scorematrix(matrix(consensus[k, ], 1, N), use_cpp = use_cpp)
      if (is.null(Wk)) {
        TauX[k] <- sum(cij * Sij) / (nj * (N * (N - 1)))
      } else {
        TauX[k] <- sum(cij * Sij) / (sum(Wk) * (N * (N - 1)))
      }
    }
  } else {
    Sij <- scorematrix(consensus, use_cpp = use_cpp)
    if (is.null(Wk)) {
      TauX <- sum(cij * Sij) / (nj * (N * (N - 1)))
    } else {
      TauX <- sum(cij * Sij) / (sum(Wk) * (N * (N - 1)))
    }
  }
  
  return(list(Consensus = reordering(consensus), Tau = TauX))
}


# ══════════════════════════════════════════════════════════════════════════
# BBconsensus3 - MAIN BB LOOP (ORIGINAL LOGIC)
# ══════════════════════════════════════════════════════════════════════════

BBconsensus3 <- function(RR, cij, Po, PS, buckets, use_cpp) {
  
  CR <- RR
  a <- t(matrix(sort(RR)))
  ord <- t(matrix(order(RR)))
  r <- ReorderingBB(RR, use_cpp = use_cpp)
  
  BR.R <- r
  BR.P <- 0
  WCi <- 1
  lambda <- 1
  N <- ncol(RR)
  
  while (WCi == 1) {
    
    if (PS) message("BB Round ", lambda)
    
    # PRIMARY LOOP: add k-th object
    for (k in 2:ncol(a)) {
      
      B <- nrow(BR.R)
      b <- 1:k
      
      # Initialize consolidation variables
      KR.R <- NULL
      KR.P <- NULL
      
      # SECONDARY LOOP: check branches
      for (nb in 1:B) {
        
        BR.R[nb, ] <- ReorderingBB(t(matrix(BR.R[nb, ])), use_cpp = use_cpp)
        
        # Generate and filter branches
        rpbr <- branchesbuckets(
          matrix(BR.R[nb, ], nrow = 1), cij, b, Po, ord,
          matrix(BR.P[nb]), buckets, use_cpp
        )
        
        R <- rpbr$cR
        Pbr <- rpbr$pcR
        
        # ═══════════════════════════════════════════════════════════
        # RATTAPPO LOGIC - ORIGINAL
        # ═══════════════════════════════════════════════════════════
        
        if (!is.null(R) && length(R) > 0) {
          
          rattappo <- matrix(0, nrow(R), 1)
          
          if (k <= buckets) {
            # Early phase: remove if > buckets
            for (L in 1:nrow(R)) {
              rattappo[L, 1] <- ifelse(length(table(R[L, ord[b]])) > buckets, 1, 0)
            }
            
          } else {
            # Late phase: remove if > buckets OR < (buckets-2)
            for (L in 1:nrow(R)) {
              rattappo[L, 1] <- ifelse(
                length(table(R[L, ord[b]])) > buckets | 
                  length(table(R[L, ord[b]])) < (buckets - 2),
                1, 0
              )
            }
          }
          
          # Apply filter
          if (sum(rattappo > 0)) {
            valid_idx <- which(rattappo == 0)
            if (length(valid_idx) > 0) {
              R <- R[valid_idx, , drop = FALSE]
              Pbr <- matrix(Pbr[valid_idx], ncol = 1)
            } else {
              R <- NULL
              Pbr <- NULL
            }
          }
          
          if (!is.null(R) && is.null(nrow(R))) {
            R <- matrix(R, nrow = 1)
          }
        }
        
        # ═══════════════════════════════════════════════════════════
        # Consolidate branches
        # ═══════════════════════════════════════════════════════════
        
        if (!is.null(R) && length(R) > 0) {
          if (is.null(KR.R)) {
            KR.R <- R
            KR.P <- Pbr
          } else {
            KR.R <- rbind(KR.R, R)
            KR.P <- rbind(KR.P, Pbr)
          }
        }
        
      }  # end secondary loop
      
      # Update branches for next k
      if (!is.null(KR.R) && length(KR.R) > 0) {
        BR.R <- KR.R
        BR.P <- matrix(KR.P, ncol = 1)
      }
      
      if (PS && !is.null(BR.R)) {
        message("  k=", k, ": ", nrow(BR.R), " branches")
      }
      
    }  # end primary loop
    
    # ═══════════════════════════════════════════════════════════════
    # CHECK AND CONVERGENCE
    # ═══════════════════════════════════════════════════════════════
    
    # Remove solutions with wrong number of buckets
    if (!is.null(BR.R) && nrow(BR.R) > 0) {
      checklength <- which(apply(BR.R, 1, function(x) length(table(x))) != buckets)
      
      if (length(checklength) > 0 && length(checklength) < nrow(BR.R)) {
        BR.R <- matrix(BR.R[-checklength, ], ncol = N)
        BR.P <- matrix(BR.P[-checklength], ncol = 1)
      } else if (length(checklength) == nrow(BR.R)) {
        # All solutions have wrong bucket count - this shouldn't happen
        warning("No valid solutions found with exactly ", buckets, " buckets")
        WCi <- 0
        next
      }
    }
    
    SSP <- matrix(which(BR.P == min(BR.P)))
    MinP <- min(BR.P)
    PenMin <- Po - MinP
    
    # Check convergence
    if (PenMin == 0) {
      # Found optimal solution
      CR <- matrix(BR.R[SSP, ], length(SSP), N)
      WCi <- 0
      
    } else {
      # Continue to next round
      Po <- MinP
      WCi <- 1
      lambda <- lambda + 1
      nRR <- matrix(BR.R[SSP[1], ], 1, N)
      BR.R <- nRR
      BR.P <- 0
      a <- t(matrix(sort(BR.R)))
      ord <- t(matrix(order(BR.R)))
    }
    
  }  # end while
  
  return(CR)
}


# ══════════════════════════════════════════════════════════════════════════
# branchesbuckets - Generate and filter branches
# ══════════════════════════════════════════════════════════════════════════

branchesbuckets <- function(brR, cij, b, Po, ord, Pb, buckets, use_cpp) {
  
  # Generate candidate branches
  candidate <- findbranchesbuckets(brR, ord, b, buckets)
  
  if (is.null(candidate) || nrow(candidate) == 0) {
    return(list(cR = NULL, pcR = NULL))
  }
  
  Pb <- matrix(rep(Pb, nrow(candidate)))
  
  # ══════════════════════════════════════════════════════════════════════════
  # OTTIMIZZAZIONE C++: Calcolo batch delle penalità
  # ══════════════════════════════════════════════════════════════════════════
  
  if (use_cpp) {
    # Use C++ batch penalty calculation
    addpenalty <- PenaltyBB2_batch_impl(cij, candidate, as.integer(ord[b]))
  } else {
    # R fallback
    addpenalty <- matrix(0, nrow(candidate), 1)
    for (gm in 1:nrow(candidate)) {
      addpenalty[gm, ] <- PenaltyBB2(cij, candidate[gm, ], ord[b], use_cpp = FALSE)
    }
  }
  
  # Filter by penalty
  Pbr <- addpenalty + Pb
  idp <- which(Pbr <= Po)
  
  if (length(idp) == 0) {
    return(list(cR = NULL, pcR = NULL))
  }
  
  R <- candidate[idp, , drop = FALSE]
  Pbr <- matrix(Pbr[idp], ncol = 1)
  
  if (is.null(nrow(R))) {
    R <- matrix(R, nrow = 1)
  }
  
  return(list(cR = R, pcR = Pbr))
}


# ══════════════════════════════════════════════════════════════════════════
# findbranchesbuckets - Generate branches (ORIGINAL LOGIC)
# ══════════════════════════════════════════════════════════════════════════

findbranchesbuckets <- function(R, ord, b, buckets) {
  
  KR <- R[ord[b]]
  KR <- KR[-length(KR)]
  mo <- max(KR)
  mi <- min(KR)
  KR[length(KR) + 1] <- mo + 1
  R[ord[b]] <- KR
  candidate <- matrix(R, nrow = 1)
  KO <- 1
  
  while (KO == 1) {
    
    R[ord[b[length(b)]]] <- R[ord[b[length(b)]]] - 1
    
    if (length(table(R[ord[b]])) <= (buckets + 1)) {
      candidate <- rbind(candidate, R)
    }
    
    if (mi - R[ord[b[length(b)]]] > 1) {
      KO <- 0
    }
  }
  
  # Remove duplicates (original logic)
  Rt <- candidate[!duplicated(reordering(candidate[, ord[b], drop = FALSE])), , drop = FALSE]
  
  return(Rt)
}


# ══════════════════════════════════════════════════════════════════════════
# BBconsensusBuckts - Phase 1 initialization (OPTIMIZED)
# ══════════════════════════════════════════════════════════════════════════

BBconsensusBuckts <- function(RR, cij, PS, buckets, use_cpp) {
  
  if (is.null(ncol(RR))) {
    RR <- matrix(RR, nrow = 1)
  }
  
  CR <- RR
  N <- ncol(RR)
  sij <- scorematrix(RR, use_cpp = use_cpp)
  Po <- sum(abs(cij)) - sum(cij * sij)
  a <- t(matrix(sort(RR, decreasing = TRUE)))
  ord <- t(matrix(order(RR, decreasing = TRUE)))
  R <- RR
  addpenalty <- matrix(0, length(a), 1)
  
  # Explore initial solution
  for (k in 2:length(a)) {
    
    b <- 1:k
    R <- ReorderingBB(R, use_cpp = use_cpp)
    
    # OPTIMIZED: Direct assignment instead of t(matrix(...))
    KR <- R[ord[b]]
    KR <- KR[-length(KR)]
    mo <- max(KR)
    mi <- min(KR)
    aa <- 1
    KO <- 1
    KR[length(KR) + 1] <- mo + 1
    R[ord[b]] <- KR
    
    candidate <- matrix(0, nrow(RR), ncol(RR))
    Pb <- matrix(0, 1, 1)
    Pc <- 1
    
    while (KO == 1) {
      
      candidate <- rbind(candidate, R)
      
      if (aa == 1) {
        candidate <- matrix(candidate[-1, ], 1, ncol(candidate))
      }
      
      Sij <- scorematrix(matrix(candidate[aa, ], 1, N), use_cpp = use_cpp)
      Pb <- rbind(Pb, sum(abs(cij)) - sum(cij * Sij))
      
      if (aa == 1) {
        Pb <- matrix(Pb[-1, ], 1, 1)
      }
      
      # Penalize wrong bucket count
      if (length(table(candidate[aa, ])) != buckets) {
        Pb[aa] <- 1e10
      }
      
      if (Pb[aa] == 0) {
        CR <- R
        Po <- 0
        Pc <- 0
        break
      }
      
      Pc <- 1
      R[ord[b[length(b)]]] <- R[ord[b[length(b)]]] - 1
      
      if (mi - R[ord[b[length(b)]]] > 1) {
        KO <- 0
      }
      
      aa <- aa + 1
    }
    
    if (Pc == 0) break
    
    minp <- min(Pb)
    posp <- which(Pb == min(Pb))
    
    if (minp <= Po) {
      Po <- minp
      CR <- t(matrix(candidate[posp[1], ]))
      R <- CR
      addpenalty[k, 1] <- PenaltyBB2(cij, R, ord[b], use_cpp = use_cpp)
    } else {
      R <- CR
      addpenalty[k, 1] <- PenaltyBB2(cij, R, ord[b], use_cpp = use_cpp)
    }
    
    candidate <- matrix(0, nrow(R), N)
    Pb <- matrix(0, 1, 1)
  }
  
  if (Pc == 0) {
    Po <- 0
    addpenalty <- 0
  }
  
  Po <- sum(addpenalty)
  
  return(list(cons = CR, pen = Po))
}


# ══════════════════════════════════════════════════════════════════════════
# N2R - Generate random initial candidates
# ══════════════════════════════════════════════════════════════════════════

N2R <- function(N, buckets, cij, M, repts, use_cpp = TRUE) {
  
  R <- matrix(0, repts, N)
  ta <- numeric(length = repts)
  
  for (j in 1:repts) {
    
    repeat {
      r <- round(runif(N, 1, buckets))
      lt <- length(table(r))
      
      if (lt == buckets) {
        t <- sum(scorematrix(r, use_cpp = use_cpp) * cij / (M * N * (N - 1)))
        if (t > 0) break
      }
    }
    
    R[j, ] <- r
    ta[j] <- t
  }
  
  # Remove duplicates
  if (sum(duplicated(R)) != 0) {
    idr <- which(duplicated(R))
    R <- R[-idr, , drop = FALSE]
    ta <- ta[-idr]
  }
  
  # Sort by tau (descending)
  id <- order(ta, decreasing = TRUE)
  R <- R[id, , drop = FALSE]
  ta <- ta[id]
  
  return(list(R = R, ta = ta))
}


# ══════════════════════════════════════════════════════════════════════════
# Quick algorithm for Median Constrained Bucket Order
# ══════════════════════════════════════════════════════════════════════════

Quick_buckets <- function(cij, nbuckets, nj, Wk = NULL, maxiter = 10, PS = TRUE, use_cpp = TRUE) {
  
  N <- ncol(cij)
  isw <- is.null(Wk)
  
  if (PS) {
    message("Quick algorithm: generating ", maxiter, " candidate solutions...")
  }
  
  # ══════════════════════════════════════════════════════════════════════════
  # STEP 1: Generate multiple random starting points
  # ══════════════════════════════════════════════════════════════════════════
  
  all_solutions <- vector("list", maxiter)
  all_taus <- numeric(maxiter)
  
  for (iter in 1:maxiter) {
    
    if (PS) {
      message("  Iteration ", iter, "/", maxiter)
    }
    
    # Generate random initial candidate with nbuckets
    R_init <- findconsensusBB_buckets(cij, nbuckets)
    
    # Optionally perturb the initial solution for diversity
    if (iter > 1) {
      # Add small random perturbation
      perm_idx <- sample(N, size = min(3, N))
      for (idx in perm_idx) {
        # Random bucket assignment (within 1:nbuckets)
        R_init[1, idx] <- sample(1:nbuckets, 1)
      }
    }
    
    # Calculate initial penalty
    Sij_init <- scorematrix(R_init, use_cpp = use_cpp)
    Po_init <- sum(abs(cij)) - sum(cij * Sij_init)
    
    # Refine with lightweight BB (1 pass only)
    bb_result <- BBconsensusBuckts_phase1(
      RR = R_init,
      cij = cij,
      nbuckets = nbuckets,
      PS = FALSE,
      use_cpp = use_cpp
    )
    
    consensus <- bb_result$cons
    
    # Ensure it's a matrix
    if (!is.matrix(consensus)) {
      consensus <- matrix(consensus, 1, N)
    }
    
    # Calculate TauX for this solution
    Sij <- scorematrix(consensus, use_cpp = use_cpp)
    
    if (!isw) {
      tau <- sum(cij * Sij) / (sum(Wk) * (N * (N - 1)))
    } else {
      tau <- sum(cij * Sij) / (nj * (N * (N - 1)))
    }
    
    all_solutions[[iter]] <- consensus
    all_taus[iter] <- tau
    
    if (PS) {
      message("    Tau: ", round(all_taus[iter], 4))
    }
  }
  
  # ══════════════════════════════════════════════════════════════════════════
  # STEP 2: Combine and select best unique solutions
  # ══════════════════════════════════════════════════════════════════════════
  
  # Combine all solutions
  all_consensus <- do.call(rbind, all_solutions)
  
  # Remove duplicates
  unique_idx <- !duplicated(all_consensus)
  unique_consensus <- all_consensus[unique_idx, , drop = FALSE]
  unique_taus <- all_taus[unique_idx]
  
  # Keep only the best solutions
  best_tau <- max(unique_taus)
  best_idx <- which(unique_taus == best_tau)
  
  final_consensus <- unique_consensus[best_idx, , drop = FALSE]
  final_tau <- unique_taus[best_idx]
  
  if (PS) {
    message("Quick completed: ", length(best_idx), " solution(s) with Tau = ", 
            round(best_tau, 4))
  }
  
  # Apply reordering before returning
  final_consensus <- reordering(final_consensus)
  
  return(list(Consensus = final_consensus, Tau = final_tau))
}


# ══════════════════════════════════════════════════════════════════════════
# HELPER FUNCTION: Quick initialization for BB
# ══════════════════════════════════════════════════════════════════════════

findconsensusBB_buckets <- function(cij, nbuckets) {
  
  N <- ncol(cij)
  
  # Edge case: nbuckets = N-1 (almost full ranking)
  # Just assign each item its own bucket except tie the two weakest
  if (nbuckets >= N - 1) {
    strength <- colSums(cij)
    ranking <- rank(-strength, ties.method = "first")
    
    # Find two items with closest strength to tie together
    strength_sorted <- sort(strength)
    min_diff_idx <- which.min(diff(strength_sorted))
    tie_value <- strength_sorted[min_diff_idx]
    
    # Assign buckets: most items get unique bucket, two weakest share one
    bucket <- ranking
    tie_items <- which(strength <= strength_sorted[min_diff_idx + 1])
    if (length(tie_items) >= 2) {
      bucket[tie_items[1:2]] <- max(bucket)
    }
    
    # Renumber to 1:nbuckets
    bucket_map <- rank(unique(sort(bucket)))
    final_ranking <- bucket_map[match(bucket, unique(sort(bucket)))]
    
    return(matrix(final_ranking, 1, N))
  }
  
  # Normal case: nbuckets << N
  strength <- colSums(cij)
  
  # Use more robust method: assign items to buckets by strength percentiles
  # This avoids issues with identical quantile values
  strength_ranks <- rank(strength, ties.method = "first")
  
  # Assign to buckets based on rank
  items_per_bucket <- N / nbuckets
  bucket <- ceiling(strength_ranks / items_per_bucket)
  
  # Ensure exactly nbuckets
  bucket <- pmin(bucket, nbuckets)
  
  # Fix if we have too few buckets (due to ties)
  actual_buckets <- length(unique(bucket))
  
  if (actual_buckets < nbuckets) {
    # Split largest buckets
    while (length(unique(bucket)) < nbuckets) {
      bucket_sizes <- table(bucket)
      largest_bucket <- as.numeric(names(which.max(bucket_sizes)))
      items_in_largest <- which(bucket == largest_bucket)
      
      if (length(items_in_largest) < 2) break  # Can't split further
      
      # Split in half
      mid <- length(items_in_largest) %/% 2
      new_bucket_id <- max(bucket) + 1
      bucket[items_in_largest[(mid+1):length(items_in_largest)]] <- new_bucket_id
    }
  } else if (actual_buckets > nbuckets) {
    # Merge smallest buckets
    while (length(unique(bucket)) > nbuckets) {
      bucket_sizes <- table(bucket)
      smallest_bucket <- as.numeric(names(which.min(bucket_sizes)))
      
      # Merge with adjacent bucket (by strength)
      bucket_means <- tapply(strength, bucket, mean)
      other_buckets <- setdiff(unique(bucket), smallest_bucket)
      
      if (length(other_buckets) == 0) break
      
      distances <- abs(bucket_means[as.character(other_buckets)] - 
                         bucket_means[as.character(smallest_bucket)])
      adjacent <- other_buckets[which.min(distances)]
      
      bucket[bucket == smallest_bucket] <- adjacent
    }
  }
  
  # Renumber buckets 1:nbuckets by mean strength
  bucket_means <- tapply(strength, bucket, mean)
  bucket_order <- rank(bucket_means)
  
  final_ranking <- bucket_order[bucket]
  
  return(matrix(final_ranking, 1, N))
}


# ══════════════════════════════════════════════════════════════════════════
# BBconsensusBuckts_phase1: Lightweight BB for Quick (1 pass only)
# ══════════════════════════════════════════════════════════════════════════

BBconsensusBuckts_phase1 <- function(RR, cij, nbuckets, PS = FALSE, use_cpp = TRUE) {
  
  if (is.null(ncol(RR))) {
    RR <- matrix(RR, nrow = 1)
  }
  
  CR <- RR
  N <- ncol(RR)
  sij <- scorematrix(RR, use_cpp = use_cpp)
  Po <- sum(abs(cij)) - sum(cij * sij)
  a <- t(matrix(sort(RR, decreasing = TRUE)))
  ord <- t(matrix(order(RR, decreasing = TRUE)))
  R <- RR
  addpenalty <- matrix(0, length(a), 1)
  
  # Single pass exploration
  for (k in 2:length(a)) {
    
    b <- 1:k
    R <- ReorderingBB(R, use_cpp = use_cpp)
    
    # OPTIMIZED: Direct assignment
    KR <- R[ord[b]]
    KR <- KR[-length(KR)]
    mo <- max(KR)
    mi <- min(KR)
    aa <- 1
    KO <- 1
    KR[length(KR) + 1] <- mo + 1
    R[ord[b]] <- KR
    
    candidate <- matrix(0, nrow(RR), ncol(RR))
    Pb <- matrix(0, 1, 1)
    Pc <- 1
    
    while (KO == 1) {
      
      candidate <- rbind(candidate, R)
      
      if (aa == 1) {
        candidate <- matrix(candidate[-1, ], 1, ncol(candidate))
      }
      
      Sij <- scorematrix(matrix(candidate[aa, ], 1, N), use_cpp = use_cpp)
      Pb <- rbind(Pb, sum(abs(cij)) - sum(cij * Sij))
      
      if (aa == 1) {
        Pb <- matrix(Pb[-1, ], 1, 1)
      }
      
      # Penalize wrong bucket count
      if (length(table(candidate[aa, ])) != nbuckets) {
        Pb[aa] <- 1e10
      }
      
      # Check if optimal found
      if (Pb[aa] == 0) {
        CR <- R
        Po <- 0
        Pc <- 0
        break
      }
      
      R[ord[b[length(b)]]] <- R[ord[b[length(b)]]] - 1
      
      if (mi - R[ord[b[length(b)]]] > 1) {
        KO <- 0
      }
      
      aa <- aa + 1
    }
    
    if (Pc == 0) {
      break
    }
    
    # Update with best found
    minp <- min(Pb)
    posp <- which(Pb == min(Pb))
    
    if (minp <= Po) {
      Po <- minp
      CR <- t(matrix(candidate[posp[1], ]))
      R <- CR
      addpenalty[k, 1] <- PenaltyBB2(cij, R, ord[b], use_cpp = use_cpp)
    } else {
      R <- CR
      addpenalty[k, 1] <- PenaltyBB2(cij, R, ord[b], use_cpp = use_cpp)
    }
    
    candidate <- matrix(0, nrow(R), N)
    Pb <- matrix(0, 1, 1)
  }
  
  if (Pc == 0) {
    Po <- 0
  } else {
    Po <- sum(addpenalty)
  }
  
  return(list(cons = CR, pen = Po))
}

# ══════════════════════════════════════════════════════════════════════════
# DECoR algorithm for Median Constrained Bucket Order
# Differential Evolution for Consensus Ranking with bucket constraints
# Based on original code by Giulio Mazzeo and Antonio D'Ambrosio
# ══════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════
# DECoR algorithm for Median Constrained Bucket Order
# Differential Evolution for Consensus Ranking with bucket constraints
# Based on original code by Giulio Mazzeo and Antonio D'Ambrosio
# ══════════════════════════════════════════════════════════════════════════

DECoR_buckets <- function(cij, nbuckets, nj, Wk = NULL, 
                          np = 10, gl = 100, ff = 0.4, cr = 0.8,
                          maxiter = 10, PS = TRUE, use_cpp = TRUE) {
  
  tic <- proc.time()[3]
  N <- ncol(cij)
  
  solutions <- matrix(0, 1, N)
  Taos <- 0
  Cs <- 0
  
  if (PS) {
    message("DECoR algorithm: running ", maxiter, " iterations...")
  }
  
  # Run multiple DECoR iterations
  for (j in 1:maxiter) {
    
    dcb <- decorbuckets(
      NP = np,
      L = gl, 
      FF = ff,
      CR = cr,
      cij = cij,
      NJ = nj,
      buckets = nbuckets,
      use_cpp = use_cpp
    )
    
    solutions <- rbind(solutions, dcb$ConsR)
    Cs <- c(Cs, dcb$bestcost)
    Taos <- c(Taos, dcb$Tau)
    
    if (PS && (j %% 5 == 0 || j == 1)) {
      message("  Iteration ", j, "/", maxiter)
    }
  }
  
  # Remove initialization row
  solutions <- solutions[-1, , drop = FALSE]
  Taos <- Taos[-1]
  Cs <- Cs[-1]
  
  # Find best solution
  index <- which.max(Taos)
  oversol <- solutions[index, , drop = FALSE]
  TauX <- Taos[index]
  Cost <- Cs[index]
  
  if (!is.matrix(oversol)) {
    oversol <- matrix(oversol, nrow = 1)
  }
  
  toc <- proc.time()[3]
  eltime <- toc - tic
  
  if (PS) {
    message("DECoR completed: Best Tau = ", round(TauX, 4))
  }
  
  return(list(
    Consensus = reordering(oversol), 
    Tau = TauX, 
    Eltime = eltime
  ))
}


# ══════════════════════════════════════════════════════════════════════════
# Single DECoR run - one evolutionary cycle
# ══════════════════════════════════════════════════════════════════════════

decorbuckets <- function(NP, L, FF, CR, cij, NJ, buckets, use_cpp = TRUE) {
  
  tic <- proc.time()[3]
  N <- nrow(cij)
  
  # ════════════════════════════════════════════════════════════════════════
  # STEP 1: Initialize population with valid bucket solutions
  # ════════════════════════════════════════════════════════════════════════
  
  # EFFICIENT population generation - no infinite loops!
  genpop <- function(N, buckets) {
    # Start with guaranteed coverage of all buckets
    x <- rep(1:buckets, length.out = N)
    # Shuffle to randomize
    x <- sample(x)
    return(x)
  }
  
  population <- t(replicate(NP, genpop(N, buckets)))
  
  # ════════════════════════════════════════════════════════════════════════
  # STEP 2: Evaluate initial population
  # ════════════════════════════════════════════════════════════════════════
  
  costs <- taos <- matrix(0, NP, 1)
  
  for (i in 1:NP) {
    COTA <- combincost(population[i, ], cij, NJ, use_cpp)
    costs[i] <- COTA$cp
    taos[i] <- COTA$tp
  }
  
  # Store best individual
  bestc <- min(costs)
  bestind <- which(costs == min(costs))
  bestT <- max(taos)
  besti <- population[bestind, , drop = FALSE]
  
  # ════════════════════════════════════════════════════════════════════════
  # STEP 3: Evolution loop
  # ════════════════════════════════════════════════════════════════════════
  
  g <- 2
  no_gain <- 0
  
  while (no_gain < L) {
    
    # Individuals mutation
    for (i in 1:NP) {
      
      # MUTATION: rand/1/bin
      evolution <- mutaterand1buckets(population, FF, i, buckets)
      
      # CROSSOVER
      evolution <- crossoverbuckets(population[i, ], evolution, CR, buckets)
      
      # REPAIR: Discretization and bucket correction
      evolution <- childclosintbuckets(evolution, buckets)
      
      # CRITICAL VALIDATION: Ensure exactly 'buckets' unique values
      if (length(unique(evolution)) != buckets) {
        # Repair failed - regenerate valid solution
        evolution <- genpop(N, buckets)
      }
      
      # EVALUATE
      COTAN <- combincost(evolution, cij, NJ, use_cpp)
      cost_new <- COTAN$cp
      ta_new <- COTAN$tp
      
      # SELECTION: Keep better solution
      if (cost_new < costs[i]) {
        population[i, ] <- evolution
        costs[i] <- cost_new
        taos[i] <- ta_new
      }
    }
    
    # Store best individual of current generation
    bestco <- min(costs)
    bestc <- c(bestc, bestco)
    bestind <- which.min(costs)
    bestTa <- max(taos)
    bestT <- c(bestT, bestTa)
    bestin <- population[bestind, ]
    besti <- rbind(besti, bestin)
    
    # Check if this generation improved solutions
    if (bestc[g] == bestc[(g - 1)]) {
      no_gain <- no_gain + 1
    } else {
      no_gain <- 0
    }
    
    g <- g + 1
  }
  
  # ════════════════════════════════════════════════════════════════════════
  # STEP 4: Select ALL best solutions (with bucket validation)
  # ════════════════════════════════════════════════════════════════════════
  
  indexes <- which(bestc == min(bestc))
  
  if (length(indexes) == 1) {
    bests <- reordering(matrix(besti[indexes, ], 1, N))
  } else {
    bests <- reordering(besti[indexes, , drop = FALSE])
  }
  
  # CRITICAL: Filter to keep only solutions with exactly 'buckets'
  valid_idx <- apply(bests, 1, function(x) length(unique(x)) == buckets)
  
  if (sum(valid_idx) == 0) {
    # No valid solutions found - regenerate one
    warning("No solutions with exactly ", buckets, " buckets found. Regenerating...")
    bests <- matrix(genpop(N, buckets), 1, N)
    avgTau <- bestT[indexes[1]]
  } else {
    bests <- bests[valid_idx, , drop = FALSE]
    avgTau <- bestT[indexes[valid_idx]]
  }
  
  ConsR <- unique(bests)
  
  mycheck <- unique(avgTau)
  if (length(mycheck) > 1) {
    warning("Multiple Tau values found in best solutions")
  }
  
  Tau <- matrix(rep(mycheck, nrow(ConsR)), nrow(ConsR), 1)
  bestcost <- unique(bestc[indexes])
  
  toc <- proc.time()[3]
  eltime <- toc - tic
  
  return(list(
    ConsR = ConsR,
    Tau = Tau,
    besti = besti,
    bestc = bestc,
    bestcost = bestcost,
    bests = bests,
    avgTau = avgTau,
    bestT = bestT,
    Eltime = eltime
  ))
}


# ══════════════════════════════════════════════════════════════════════════
# MUTATION: rand/1/bin (ORIGINAL LOGIC)
# ══════════════════════════════════════════════════════════════════════════

mutaterand1buckets <- function(X, FF, i, buckets) {
  
  D <- nrow(X)
  
  # Select 3 random individuals different from current
  a <- sample(D)
  
  for (j in 1:3) {
    if (a[j] == i) {
      a[j] <- a[4]
    }
  }
  
  r1 <- a[1]
  r2 <- a[2]
  r3 <- a[3]
  
  # Apply mutation (Rand1Bin)
  v <- X[r1, ] + FF * (X[r2, ] - X[r3, ])
  
  # Bound to valid range
  v[v >= (buckets + 0.5)] <- buckets + 0.4
  v[v <= 0.4] <- 0.5
  
  return(v)
}


# ══════════════════════════════════════════════════════════════════════════
# CROSSOVER (ORIGINAL LOGIC)
# ══════════════════════════════════════════════════════════════════════════

crossoverbuckets <- function(x, v, CR, buckets) {
  
  if (!is.matrix(x)) {
    x <- matrix(x, 1, length(x))
  }
  
  D <- ncol(x)
  u <- matrix(0, 1, D)
  
  for (i in 1:D) {
    if (runif(1) > CR) {
      u[i] <- v[i]
    } else {
      u[i] <- x[i]
    }
  }
  
  # Bound to valid range
  u[u > buckets] <- buckets
  u[u < 1] <- 1
  
  return(u)
}


# ══════════════════════════════════════════════════════════════════════════
# REPAIR: Ensure exactly nbuckets unique values (ORIGINAL LOGIC)
# ══════════════════════════════════════════════════════════════════════════

childclosintbuckets <- function(r, buckets) {
  
  D <- length(r)
  
  # Closest integer
  x <- round(r)
  
  # Correct out of bound
  out_of_bounds <- which(x > buckets | x < 1)
  if (length(out_of_bounds) > 0) {
    x[out_of_bounds] <- sample(1:buckets, length(out_of_bounds), replace = TRUE)
  }
  
  # ════════════════════════════════════════════════════════════════════════
  # CRITICAL: Ensure exactly 'buckets' unique values - EFFICIENT VERSION
  # ════════════════════════════════════════════════════════════════════════
  
  # Find missing buckets
  present_buckets <- unique(x)
  missing_buckets <- setdiff(1:buckets, present_buckets)
  
  if (length(missing_buckets) == 0) {
    # Already has exactly 'buckets' unique values
    return(x)
  }
  
  # Find duplicates to replace
  dup_idx <- which(duplicated(x))
  
  if (length(dup_idx) == 0) {
    # No duplicates but missing buckets - force random positions
    positions_to_replace <- sample(D, length(missing_buckets))
    x[positions_to_replace] <- missing_buckets
    return(x)
  }
  
  # ════════════════════════════════════════════════════════════════════════
  # EFFICIENT ASSIGNMENT - No loops!
  # ════════════════════════════════════════════════════════════════════════
  
  n_missing <- length(missing_buckets)
  n_dup <- length(dup_idx)
  
  if (n_missing <= n_dup) {
    # Simple case: enough duplicates for all missing buckets
    x[dup_idx[1:n_missing]] <- missing_buckets
  } else {
    # Not enough duplicates: need to replace more positions
    # First, use all duplicates
    x[dup_idx] <- missing_buckets[1:n_dup]
    
    # Then force remaining missing buckets into random positions
    still_missing <- missing_buckets[(n_dup + 1):n_missing]
    available_positions <- setdiff(1:D, dup_idx)
    positions_to_replace <- sample(available_positions, length(still_missing))
    x[positions_to_replace] <- still_missing
  }
  
  return(x)
}


# ══════════════════════════════════════════════════════════════════════════
# HELPER: Compute cost and tau (ORIGINAL LOGIC with C++ optimization)
# ══════════════════════════════════════════════════════════════════════════

combincost <- function(ranking, cij, M, use_cpp = TRUE) {
  
  if (!is.matrix(ranking)) {
    ranking <- matrix(ranking, 1, length(ranking))
  }
  
  N <- ncol(ranking)
  
  # Compute score matrix (with C++ if available)
  sij <- scorematrix(ranking, use_cpp = use_cpp)
  
  # Max distance
  maxdist <- (N * (N - 1))
  
  # Computing the tau
  t <- sum(sum(cij * sij)) / (M * (N * (N - 1)))
  
  # Computing the distance (from tau)
  c <- (maxdist / 2) * (1 - t) * M
  
  return(list(tp = t, cp = c))
}