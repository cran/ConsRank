#include <Rcpp.h>
using namespace Rcpp;

 //' @noRd
 // [[Rcpp::export]]
 NumericVector PenaltyBB2_batch_impl(NumericMatrix cij, NumericMatrix candidates, IntegerVector ord) {
   
   int n_candidates = candidates.nrow();
   int ord_len = ord.size();
   NumericVector result(n_candidates);
   
   // Pre-compute sign matrix of cij for efficiency
   int n = cij.nrow();
   IntegerMatrix scij(n, n);
   for (int i = 0; i < n; i++) {
     for (int j = 0; j < n; j++) {
       if (cij(i, j) > 0) scij(i, j) = 1;
       else if (cij(i, j) < 0) scij(i, j) = -1;
       else scij(i, j) = 0;
     }
   }
   
   // Process each candidate
   for (int c = 0; c < n_candidates; c++) {
     NumericVector candidate = candidates(c, _);
     double total_penalty = 0.0;
     
     int last_ord_idx = ord[ord_len - 1] - 1;  // Convert to 0-based
     
     // Loop over ord elements (except last)
     for (int k = 0; k < ord_len - 1; k++) {
       int ord_k_idx = ord[k] - 1;  // Convert to 0-based
       
       // Compute Ds = sign(candidate[ord[last]] - candidate[ord[k]])
       double diff = candidate[last_ord_idx] - candidate[ord_k_idx];
       int Ds = (diff > 0) ? 1 : ((diff < 0) ? -1 : 0);
       
       double penalty = 0.0;
       int s_last_k = scij(last_ord_idx, ord_k_idx);
       int s_k_last = scij(ord_k_idx, last_ord_idx);
       
       // Complex nested if-else logic from PenaltyBB2
       if (Ds == 1) {
         if (s_last_k == 1 && s_k_last == -1) {
           penalty = cij(last_ord_idx, ord_k_idx) - cij(ord_k_idx, last_ord_idx);
         } else if ((s_last_k == 1 && s_k_last == 1) ||
           (s_last_k == 0 && s_k_last == 0) ||
           (s_last_k == 1 && s_k_last == 0) ||
           (s_last_k == 0 && s_k_last == 1)) {
           penalty = cij(last_ord_idx, ord_k_idx);
         } else if (s_last_k == -1 && s_k_last == 1) {
           penalty = 0.0;
         }
       }
       else if (Ds == -1) {
         if (s_last_k == 1 && s_k_last == -1) {
           penalty = 0.0;
         } else if ((s_last_k == 1 && s_k_last == 1) ||
           (s_last_k == 0 && s_k_last == 0) ||
           (s_last_k == 1 && s_k_last == 0) ||
           (s_last_k == 0 && s_k_last == 1)) {
           penalty = cij(ord_k_idx, last_ord_idx);
         } else if (s_last_k == -1 && s_k_last == 1) {
           penalty = cij(ord_k_idx, last_ord_idx) - cij(last_ord_idx, ord_k_idx);
         }
       }
       else if (Ds == 0) {
         if (s_last_k == 1 && s_k_last == -1) {
           penalty = -cij(ord_k_idx, last_ord_idx);
         } else if ((s_last_k == 1 && s_k_last == 1) ||
           (s_last_k == 0 && s_k_last == 0) ||
           (s_last_k == 1 && s_k_last == 0) ||
           (s_last_k == 0 && s_k_last == 1)) {
           penalty = 0.0;
         } else if (s_last_k == -1 && s_k_last == 1) {
           penalty = -cij(last_ord_idx, ord_k_idx);
         }
       }
       
       total_penalty += penalty;
     }
     
     result[c] = total_penalty;
   }
   
   return result;
 }
