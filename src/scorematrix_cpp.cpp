#include <Rcpp.h>
using namespace Rcpp;

 //' @noRd
 // [[Rcpp::export]]
 NumericMatrix scorematrix_cpp(SEXP X) {
   NumericVector x;
   
   // Handle both vector and matrix input
   if (Rf_isMatrix(X)) {
     NumericMatrix mat(X);
     if (mat.nrow() != 1) {
       stop("Matrix input must have exactly 1 row");
     }
     x = mat(0, _);
   } else {
     x = as<NumericVector>(X);
   }
   
   int m = x.size();
   NumericMatrix sm(m, m);
   
   // Initialize with -1
   std::fill(sm.begin(), sm.end(), -1.0);
   
   // Compute pairwise comparisons
   for (int i = 0; i < m; i++) {
     // Diagonal
     sm(i, i) = 0.0;
     
     // Skip if NA
     if (ISNAN(x[i])) {
       for (int j = 0; j < m; j++) {
         sm(i, j) = 0.0;
       }
       continue;
     }
     
     // Compare with other elements
     for (int j = 0; j < m; j++) {
       if (i == j) continue;
       
       if (!ISNAN(x[j])) {
         sm(i, j) = (x[i] <= x[j]) ? 1.0 : -1.0;
       } else {
         sm(i, j) = 0.0;
       }
     }
   }
   
   return sm;
 }
