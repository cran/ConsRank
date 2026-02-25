#include <Rcpp.h>
using namespace Rcpp;

 //' @noRd
 // [[Rcpp::export]]
 double PenaltyBB2_impl(NumericMatrix cij, NumericVector candidate, IntegerVector ord) {
   int n = ord.size();
   double total_penalty = 0.0;
   
   // Get the last element index (C++ is 0-indexed, R is 1-indexed)
   int last_idx = ord[n - 1] - 1;  // Convert R index to C++ index
   double last_rank = candidate[last_idx];
   
   // Pre-compute sign matrix for efficiency
   NumericMatrix sign_cij(cij.nrow(), cij.ncol());
   for (int i = 0; i < cij.nrow(); i++) {
     for (int j = 0; j < cij.ncol(); j++) {
       if (cij(i, j) > 0) sign_cij(i, j) = 1;
       else if (cij(i, j) < 0) sign_cij(i, j) = -1;
       else sign_cij(i, j) = 0;
     }
   }
   
   // Loop through all elements except the last
   for (int k = 0; k < n - 1; k++) {
     int curr_idx = ord[k] - 1;  // Convert R index to C++ index
     double curr_rank = candidate[curr_idx];
     
     // Compute Ds: sign of difference
     double Ds = 0.0;
     if (last_rank > curr_rank) Ds = 1.0;
     else if (last_rank < curr_rank) Ds = -1.0;
     else Ds = 0.0;
     
     // Get relevant cij values
     double c_last_curr = cij(last_idx, curr_idx);
     double c_curr_last = cij(curr_idx, last_idx);
     int sign_last_curr = sign_cij(last_idx, curr_idx);
     int sign_curr_last = sign_cij(curr_idx, last_idx);
     
     double penalty = 0.0;
     
     if (Ds == 1.0) {
       // Case: last_rank > curr_rank
       if (sign_last_curr == 1 && sign_curr_last == -1) {
         penalty = c_last_curr - c_curr_last;
       } else if ((sign_last_curr == 1 && sign_curr_last == 1) ||
         (sign_last_curr == 0 && sign_curr_last == 0) ||
         (sign_last_curr == 1 && sign_curr_last == 0) ||
         (sign_last_curr == 0 && sign_curr_last == 1)) {
         penalty = c_last_curr;
       } else if (sign_last_curr == -1 && sign_curr_last == 1) {
         penalty = 0.0;
       }
     } else if (Ds == -1.0) {
       // Case: last_rank < curr_rank
       if (sign_last_curr == 1 && sign_curr_last == -1) {
         penalty = 0.0;
       } else if ((sign_last_curr == 1 && sign_curr_last == 1) ||
         (sign_last_curr == 0 && sign_curr_last == 0) ||
         (sign_last_curr == 1 && sign_curr_last == 0) ||
         (sign_last_curr == 0 && sign_curr_last == 1)) {
         penalty = c_curr_last;
       } else if (sign_last_curr == -1 && sign_curr_last == 1) {
         penalty = c_curr_last - c_last_curr;
       }
     } else {
       // Case: last_rank == curr_rank (tie)
       if (sign_last_curr == 1 && sign_curr_last == -1) {
         penalty = -c_curr_last;
       } else if ((sign_last_curr == 1 && sign_curr_last == 1) ||
         (sign_last_curr == 0 && sign_curr_last == 0) ||
         (sign_last_curr == 1 && sign_curr_last == 0) ||
         (sign_last_curr == 0 && sign_curr_last == 1)) {
         penalty = 0.0;
       } else if (sign_last_curr == -1 && sign_curr_last == 1) {
         penalty = -c_last_curr;
       }
     }
     
     total_penalty += penalty;
   }
   
   return total_penalty;
 }
 
 
 //' Score matrix calculation (C++ implementation)
 //' 
 //' @param x Ranking vector
 //' @return M x M score matrix
 //' @keywords internal
 // [[Rcpp::export]]
 NumericMatrix scorematrix_impl(NumericVector x) {
   int m = x.size();
   NumericMatrix sm(m, m);
   
   // Initialize with -1
   for (int i = 0; i < m; i++) {
     for (int j = 0; j < m; j++) {
       sm(i, j) = -1.0;
     }
   }
   
   // Compute comparisons
   for (int i = 0; i < m; i++) {
     for (int j = 0; j < m; j++) {
       if (i == j) {
         sm(i, j) = 0.0;  // Diagonal
       } else if (!ISNAN(x[i]) && !ISNAN(x[j])) {
         if (x[i] <= x[j]) {
           sm(i, j) = 1.0;
         }
       } else {
         sm(i, j) = 0.0;  // NA handling
       }
     }
   }
   
   return sm;
 }
