#include "shared_functions.h"
using namespace Rcpp;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the classical GFCM
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the generalized membership matrix
//' @name calcFGCMBelongMatrix
//' @description Calculate the generalized membership matrix according to a set of
//' centroids, the observed data, the fuzziness degree, and a beta parameter
//'
//' @param centers A matrix representing the centers of the
//'   clusters with p columns and k rows
//' @param data A matrix representing the observed data with n rows
//'   and p columns
//' @param m A float representing the fuzziness degree
//' @param beta A float for the beta parameter (control speed convergence and classification crispness)
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcFGCMBelongMatrix(NumericMatrix centers, NumericMatrix data, double m, double beta){

  //calculating euclidean distance between each observation and each center
  NumericMatrix centerdistances(data.nrow(), centers.nrow());
  int i;
  for(i = 0 ; i < centers.nrow() ; i++){
    centerdistances(_,i) = calcEuclideanDistance2(data, centers(i,_));
  }

  //calculating aj
  NumericVector aj = beta * rowmins_mat(centerdistances);

  int k = centers.nrow();
  int j;
  double p = (1.0/(m-1.0));
  NumericMatrix uij(data.nrow(),centers.nrow());
  NumericVector denom(data.nrow());

  for(i = 0 ; i < k ; i++){
    NumericVector dij_aj = centerdistances(_,i) - aj;
    for(j = 0 ; j < k ; j++){
      denom = denom + pow((dij_aj) / (centerdistances(_,j) - aj),p);
    }
    NumericVector vals = (1.0/denom);
    vals[is_na(vals)] = 1;
    uij(_,i) = vals;
    denom.fill(0);
  }

  return uij;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the spatial GFCM
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the generalized membership matrix (spatial version)
//' @name calcSFGCMBelongMatrix
//' @description Calculate the generalized membership matrix (spatial version)
//' @param centers A matrix representing the centers of the clusters with p
//'   columns and k rows
//' @param data A matrix representing the observed data with n rows and p columns
//' @param wdata A matrix representing the lagged observed data with n rows
//'   and p columns
//' @param m A float representing the fuzziness degree
//' @param alpha A float representing the weight of the space in the analysis (0
//'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
//'   dimensions, 2 is twice the weight for space)
//' @param beta A float for the beta parameter (control speed convergence and classification crispness)
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcSFGCMBelongMatrix(NumericMatrix centers, NumericMatrix data, NumericMatrix wdata, double m, double alpha, double beta){

  //calculating euclidean distance between each observation and each center
  //and for the spatially lagged matrix
  NumericMatrix centerdistances(data.nrow(), centers.nrow());
  NumericMatrix Wcenterdistances(data.nrow(), centers.nrow());
  int i;
  for(i = 0 ; i < centers.nrow() ; i++){
    centerdistances(_,i) = calcEuclideanDistance2(data, centers(i,_));
    Wcenterdistances(_,i) = calcEuclideanDistance2(wdata, centers(i,_));
  }

  //calculating aj
  NumericVector aj = beta * rowmins_mat(centerdistances);

  int k = centers.nrow();
  int j;
  double p = (1.0/(m-1.0));
  NumericMatrix uij(data.nrow(),centers.nrow());
  NumericVector denom(data.nrow());

  for(i = 0 ; i < k ; i++){
    NumericVector dij_aj = centerdistances(_,i) - aj;
    NumericVector wdij = Wcenterdistances(_,i);
    for(j = 0 ; j < k ; j++){
      denom = denom + pow(
        (dij_aj + alpha * wdij) /
          (centerdistances(_,j) - aj + alpha * Wcenterdistances(_,j))
        ,p);
    }
    NumericVector vals = (1.0/denom);
    vals[is_na(vals)] = 1;
    uij(_,i) = vals;
    denom.fill(0);
  }

  return uij;
}

