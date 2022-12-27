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
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcFGCMBelongMatrix(NumericMatrix centers, NumericMatrix data, double m, double beta, NumericVector sigmas){

  //calculating euclidean distance between each observation and each center
  NumericMatrix centerdistances(data.nrow(), centers.nrow());
  int i;
  for(i = 0 ; i < centers.nrow() ; i++){
    centerdistances(_,i) = calcEuclideanDistance2(data, centers(i,_)) / sigmas(i);
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
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @param wsigmas Same as sigmas, but calculated on the spatially lagged dataset
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcSFGCMBelongMatrix(NumericMatrix centers, NumericMatrix data, NumericMatrix wdata, double m, double alpha, double beta, NumericVector sigmas, NumericVector wsigmas){

  //calculating euclidean distance between each observation and each center
  //and for the spatially lagged matrix
  NumericMatrix centerdistances(data.nrow(), centers.nrow());
  NumericMatrix Wcenterdistances(data.nrow(), centers.nrow());
  int i;
  for(i = 0 ; i < centers.nrow() ; i++){
    centerdistances(_,i) = calcEuclideanDistance2(data, centers(i,_)) / sigmas(i);
    Wcenterdistances(_,i) = calcEuclideanDistance2(wdata, centers(i,_)) /wsigmas(i);
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



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the NOISY GFCM
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the generalized membership matrix with a noise cluster
//' @name calcFGCMBelongMatrixNoisy
//' @description Calculate the generalized membership matrix according to a set of
//' centroids, the observed data, the fuzziness degree, and a beta parameter
//'
//' @param centers A matrix representing the centers of the
//'   clusters with p columns and k rows
//' @param data A matrix representing the observed data with n rows
//'   and p columns
//' @param m A float representing the fuzziness degree
//' @param beta A float for the beta parameter (control speed convergence and classification crispness)
//' @param delta A float, the value set for delta by the user
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcFGCMBelongMatrixNoisy(NumericMatrix centers, NumericMatrix data, double m, double beta, double delta, NumericVector sigmas){

  //calculating euclidean distance between each observation and each center
  int i;
  int nc = centers.nrow();
  NumericMatrix centerdistances(data.nrow(), nc+1);

  for(i = 0 ; i < centers.nrow() ; i++){
    centerdistances(_,i) = calcEuclideanDistance2(data, centers(i,_)) / sigmas(i);
  }
  NumericVector delta2(data.nrow(), delta);
  centerdistances(_,nc) = delta2;

  //calculating aj
  NumericVector aj = beta * rowmins_mat(centerdistances);

  int j;
  double p = (1.0/(m-1.0));
  NumericMatrix uij(data.nrow(), nc+1);
  NumericVector denom(data.nrow());

  for(i = 0 ; i < nc+1 ; i++){
    NumericVector dij_aj = centerdistances(_,i) - aj;
    for(j = 0 ; j < nc+1 ; j++){
      denom = denom + pow((dij_aj) / (centerdistances(_,j) - aj),p);
    }
    NumericVector vals = (1.0/denom);
    vals[is_na(vals)] = 1;
    uij(_,i) = vals;
    denom.fill(0);
  }

  NumericMatrix U2 = uij(_,Range(0,nc-1));
  return U2;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the spatial and NOISY GFCM
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the generalized membership matrix (spatial version) with a noise cluster
//' @name calcSFGCMBelongMatrixNoisy
//' @description Calculate the generalized membership matrix (spatial version) with a noise cluster
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
//' @param delta A float, the value set for delta by the user
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @param wsigmas Same as sigmas, but calculated on the spatially lagged dataset
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcSFGCMBelongMatrixNoisy(NumericMatrix centers, NumericMatrix data, NumericMatrix wdata, double m, double alpha, double beta, double delta, NumericVector sigmas, NumericVector wsigmas){

  //calculating euclidean distance between each observation and each center
  //and for the spatially lagged matrix

  int nc = centers.nrow();

  NumericMatrix centerdistances(data.nrow(), nc +1);
  NumericMatrix Wcenterdistances(data.nrow(), nc +1);
  int i;
  for(i = 0 ; i < nc ; i++){
    centerdistances(_,i) = calcEuclideanDistance2(data, centers(i,_)) / sigmas(i);
    Wcenterdistances(_,i) = calcEuclideanDistance2(wdata, centers(i,_)) /wsigmas(i);
  }

  NumericVector delta2(data.nrow(), delta);
  centerdistances(_,nc) = delta2;
  Wcenterdistances(_,nc) = delta2;

  //calculating aj
  NumericVector aj = beta * rowmins_mat(centerdistances);

  int j;
  double p = (1.0/(m-1.0));
  NumericMatrix uij(data.nrow(), nc+1);
  NumericVector denom(data.nrow());

  for(i = 0 ; i < nc+1 ; i++){
    NumericVector dij_aj = centerdistances(_,i) - aj;
    NumericVector wdij = Wcenterdistances(_,i);
    for(j = 0 ; j < nc+1 ; j++){
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

  NumericMatrix U2 = uij(_,Range(0,nc-1));
  return U2;
}

