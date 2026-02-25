#include <Rcpp.h>
using namespace Rcpp;

//' @noRd
 // [[Rcpp::export]]
 NumericMatrix kemenydesign_cpp(NumericMatrix X) {
   int n = X.nrow();
   int m = X.ncol();
   int ncols = (m * (m - 1)) / 2;
   
   NumericMatrix KX(n, ncols);
   
   for (int row = 0; row < n; row++) {
     int col_idx = 0;
     
     for (int i = 0; i < m - 1; i++) {
       for (int j = i + 1; j < m; j++) {
         double diff = X(row, i) - X(row, j);
         
         if (diff < 0) {
           KX(row, col_idx) = 1.0;
         } else if (diff > 0) {
           KX(row, col_idx) = -1.0;
         } else {
           KX(row, col_idx) = 0.0;
         }
         
         col_idx++;
       }
     }
   }
   
   return KX;
 }
