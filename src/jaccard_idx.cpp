#include "shared_functions.h"
#include <cmath>
using namespace Rcpp;



//' @title Jaccard similarity coefficient
//' @name calc_jaccard_idx
//' @description Calculate the Jaccard similarity coefficient
//' @param x A vector of positive reals
//' @param y A vector of positive reals
//' @return A double: the Jaccard similarity coefficient
//' @keywords internal
// [[Rcpp::export]]
double calc_jaccard_idx(arma::vec x, arma::vec y){
  return arma::accu(arma::min(x,y))/arma::accu(arma::max(x,y));
}


//' @title Jaccard similarity coefficient between columns of two matrices
//' @name calc_jaccard_mat
//' @description Calculate the Jaccard similarity coefficient between the
//' columns of two matrices
//' @param matX A matrix
//' @param matY A matrix
//' @return A matrix with the Jaccard index values
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix calc_jaccard_mat(NumericMatrix matX, NumericMatrix matY){
  int k = matY.ncol();
  int i,j;
  // creating an empty matrix at first
  NumericMatrix jacc_mat(k,k);
  jacc_mat.fill(-1);
  for(i = 0; i < k; i++){
    for(j = 0; j < k; j++){
      if(jacc_mat(i,j)==-1){
          jacc_mat(i,j) = calc_jaccard_idx(matX(_,j),matY(_,i));
        };
    }
  }
  return jacc_mat;
}



