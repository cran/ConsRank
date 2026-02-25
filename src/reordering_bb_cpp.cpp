#include <Rcpp.h>
using namespace Rcpp;

 //' @noRd
 // [[Rcpp::export]]
 NumericVector ReorderingBB_impl(NumericVector RR) {
   int k = RR.size();
   
   // Add 1 to all elements
   NumericVector R = clone(RR);
   for (int i = 0; i < k; i++) {
     R[i] += 1.0;
   }
   
   // Get ordering (sort indices)
   IntegerVector neword = seq(0, k - 1);
   std::sort(neword.begin(), neword.end(), 
             [&R](int i, int j) { return R[i] < R[j]; });
   
   // Calculate indexing (differences)
   NumericVector indexing(k - 1);
   for (int j = 0; j < k - 1; j++) {
     indexing[j] = R[neword[j + 1]] - R[neword[j]];
   }
   
   // Check if any zeros in indexing
   bool has_zeros = false;
   for (int j = 0; j < k - 1; j++) {
     if (indexing[j] == 0.0) {
       has_zeros = true;
       break;
     }
   }
   
   if (has_zeros) {
     int J = 0;
     while (J < k - 1) {
       if (indexing[J] == 0.0) {
         R[neword[J + 1]] = R[neword[J]];
       } else if (indexing[J] > 0.0) {
         R[neword[J + 1]] = R[neword[J]] + 2.0;
       }
       J++;
     }
   } else {
     int J = 0;
     while (J < k - 1) {
       R[neword[J + 1]] = R[neword[J]] + 2.0;
       J++;
     }
   }
   
   return R;
 }
