#include <Rcpp.h>
using namespace Rcpp;

//' @noRd
// [[Rcpp::export]]
NumericMatrix findbranches_impl(NumericMatrix R, IntegerVector ord, IntegerVector b, bool FULL) {
  
  // R must be a 1-row matrix
  if (R.nrow() != 1) {
    stop("R must be a 1-row matrix");
  }
  
  int n = R.ncol();
  int b_len = b.size();
  
  // Clone R to work with
  NumericVector R_vec = R(0, _);  // Extract first row as vector
  
  // Extract KR = R[ord[b]]
  NumericVector KR(b_len);
  for (int i = 0; i < b_len; i++) {
    KR[i] = R_vec[ord[i] - 1];  // ord is 1-based, convert to 0-based
  }
  
  // Remove last element: KR <- KR[-length(KR)]
  NumericVector KR_reduced(b_len - 1);
  for (int i = 0; i < b_len - 1; i++) {
    KR_reduced[i] = KR[i];
  }
  
  double MO = Rcpp::max(KR_reduced);
  double MI = Rcpp::min(KR_reduced);
  
  // KR[length(KR)+1] <- MO+1
  NumericVector KR_extended(b_len);
  for (int i = 0; i < b_len - 1; i++) {
    KR_extended[i] = KR_reduced[i];
  }
  KR_extended[b_len - 1] = MO + 1;
  
  // R[ord[b]] <- KR_extended
  NumericVector R_work = clone(R_vec);
  for (int i = 0; i < b_len; i++) {
    R_work[ord[i] - 1] = KR_extended[i];
  }
  
  // Initialize candidate matrix
  std::vector<NumericVector> candidates;
  
  int aa = 1;
  int KO = 1;
  int step = FULL ? 2 : 1;
  
  // Main loop: while (KO==1)
  while (KO == 1) {
    
    // candidate <- rbind(candidate, R)
    candidates.push_back(clone(R_work));
    
    // if (aa==1) candidate <- matrix(candidate[-1,], 1, ncol(candidate))
    // This removes the first empty row in R, but we start with actual data
    // so we can skip this in C++
    
    // R[ord[b[length(b)]]] <- R[ord[b[length(b)]]] - step
    int last_b_idx = ord[b_len - 1] - 1;  // Convert to 0-based
    R_work[last_b_idx] = R_work[last_b_idx] - step;
    
    // if (MI - R[ord[b[length(b)]]] > 1) KO <- 0
    if (MI - R_work[last_b_idx] > 1) {
      KO = 0;
    }
    
    aa++;
  }
  
  // Convert vector of vectors to matrix
  int n_candidates = candidates.size();
  NumericMatrix result(n_candidates, n);
  
  for (int i = 0; i < n_candidates; i++) {
    for (int j = 0; j < n; j++) {
      result(i, j) = candidates[i][j];
    }
  }
  
  return result;
}
