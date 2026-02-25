#include <Rcpp.h>
using namespace Rcpp;

//' @noRd
// [[Rcpp::export]]
NumericMatrix combinpmatr_impl(NumericMatrix X, Nullable<NumericVector> Wk = R_NilValue) {
  int n = X.nrow();  // number of judges
  int m = X.ncol();  // number of objects
  
  // Initialize combined input matrix
  NumericMatrix CI(m, m);
  
  // Get weights if provided
  bool has_weights = Wk.isNotNull();
  NumericVector weights;
  if (has_weights) {
    weights = Wk.get();
  }
  
  // For each judge (row)
  for (int idx = 0; idx < n; idx++) {
    double weight = has_weights ? weights[idx] : 1.0;
    
    // For each pair of objects (i, j)
    for (int i = 0; i < m; i++) {
      double rank_i = X(idx, i);
      
      // Skip if rank_i is NA
      if (ISNAN(rank_i)) continue;
      
      for (int j = 0; j < m; j++) {
        if (i == j) continue;  // diagonal stays 0
        
        double rank_j = X(idx, j);
        
        // Skip if rank_j is NA
        if (ISNAN(rank_j)) continue;
        
        // Compute score matrix element and add to CI
        // If rank_i <= rank_j, add weight, else subtract weight
        if (rank_i <= rank_j) {
          CI(i, j) += weight;
        } else {
          CI(i, j) -= weight;
        }
      }
    }
  }
  
  // Set row and column names if present
  if (X.hasAttribute("dimnames")) {
    List dimnames = X.attr("dimnames");
    if (dimnames.size() >= 2 && !Rf_isNull(dimnames[1])) {
      CharacterVector colnames = dimnames[1];
      CI.attr("dimnames") = List::create(colnames, colnames);
    }
  }
  
  return CI;
}
