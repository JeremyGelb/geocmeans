#include "shared_functions.h"
using namespace Rcpp;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the classical FCM
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the centroids
//' @name calcCentroids
//' @description Calculate the new centroids of the clusters based on the membership matrix
//' for a classical FCM.
//' @param data A Numeric matrix representing the observed data with n rows
//'   and p columns
//' @param belongmatrix A n X k matrix giving for each observation n, its
//'   probability to belong to the cluster k
//' @param m A float representing the fuzziness degree
//' @return A a matrix with the centers calculated for each cluster
//' @keywords internal
//
// [[Rcpp::export]]
NumericMatrix calcCentroids(arma::mat data, arma::mat belongmatrix, double m){
  arma::mat powered = pow(belongmatrix, m);
  NumericMatrix centers(belongmatrix.n_cols, data.n_cols);
  int ncB = belongmatrix.n_cols;
  int ncD = data.n_cols;
  arma::vec pi;
  NumericVector center(ncD);
  double spi;
  int i;
  int j;
  for( i = 0; i < ncB ; i++){
    pi = powered.col(i);
    spi = sum(pi);
    for(j = 0 ; j < ncD ; j++){
      center(j) = (sum((data.col(j) % pi) / spi ));
    }
    centers(i,_) = center;
  }
  return centers;
}



// data(LyonIris)
// AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
// "TxChom1564","Pct_brevet","NivVieMed")
// dataset <- LyonIris@data[AnalysisFields]
// queen <- spdep::poly2nb(LyonIris,queen=TRUE)
// Wqueen <- spdep::nb2listw(queen,style="W")
// result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
// mat <- result$Belongings[1:5,]
// df <- result$Data[1:5,]

//' @title Calculate the membership matrix
//' @name calcBelongMatrix
//' @description Calculate the membership matrix according to a set of centroids, the observed
//'  data and the fuzziness degree
//' @param centers A matrix or a dataframe representing the centers of the
//'    clusters with p columns and k rows
//' @param data A dataframe or matrix representing the observed data with n rows
//'    and p columns
//' @param m A float representing the fuzziness degree
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @return A n * k matrix representing the probability of belonging of each
//'    observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcBelongMatrix(NumericMatrix centers, NumericMatrix data, double m, NumericVector sigmas){

  //calculating euclidean distance between each observation and each center
  NumericMatrix centerdistances(data.nrow(), centers.nrow());
  int i;
  for(i = 0 ; i < centers.nrow() ; i++){
    centerdistances(_,i) = calcEuclideanDistance2(data, centers(i,_)) / sigmas(i);
  }
  NumericMatrix numerator = power_mat(centerdistances, -1.0/(m-1.0));
  NumericVector denom = rowSums(numerator);
  // finally the total matrix
  NumericMatrix belongmat(numerator.rows(),numerator.cols());
  for(i = 0; i < numerator.cols(); i++){
    NumericVector vals = numerator(_,i) / denom;
    vals[is_na(vals)] = 1;
    belongmat(_,i) = vals;
  }
  return belongmat;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the spatial version of the classical FCM (SFCM)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the centroids of SFCM
//' @name calcSWFCCentroids
//' @description Calculate the new centroids of the clusters based on the membership matrix for SFCM
//' @param data A matrix representing the observed data with n rows and p columns
//' @param wdata A matrix representing the lagged observed data with nrows and p columns
//' @param belongmatrix A n X k matrix giving for each observation n, its
//'   probability to belong to the cluster k
//' @param m An integer representing the fuzziness degree
//' @param alpha A float representing the weight of the space in the analysis (0
//'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
//'   dimensions, 2 is twice the weight for space)
//' @return A n X k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcSWFCCentroids(arma::mat data, arma::mat wdata, arma::mat belongmatrix, double m, double alpha){
  arma::mat powered = arma::pow(belongmatrix, m);
  NumericMatrix centers(belongmatrix.n_cols, data.n_cols);
  arma::mat wdata_alpha = wdata * alpha;
  arma::vec pi;
  double spi;
  int ncB = belongmatrix.n_cols;
  int ncD = data.n_cols;
  int i;
  int j;
  for( i = 0; i < ncB ; i++){
    pi = powered.col(i);
    spi = accu(pi);
    NumericVector center(ncD);
    for(j = 0 ; j < ncD ; j++){
      center(j) = sum((data.col(j) + wdata_alpha.col(j)) % pi) / ((1.0+alpha) * spi);
    }
    centers(i,_) = center;
  }
  return centers;
}


//' @title Calculate the membership matrix (spatial version)
//' @name calcSFCMBelongMatrix
//' @description Calculate the membership matrix (spatial version) according to a set of
//' centroids, the observed data, the fuzziness degree a neighbouring matrix and
//' a spatial weighting term
//' @param centers A matrix or a dataframe representing the centers of the
//'   clusters with p columns and k rows
//' @param data A matrix representing the observed data with n rows and p columns
//' @param wdata A matrix representing the lagged observed data with n rows and p columns
//' @param m A float representing the fuzziness degree
//' @param alpha A float representing the weight of the space in the analysis (0
//'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
//'   dimensions, 2 is twice the weight for space)
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @param wsigmas Same as sigmas, but calculated on the spatially lagged dataset
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcSFCMBelongMatrix(NumericMatrix centers, NumericMatrix data, NumericMatrix wdata, double m, double alpha, NumericVector sigmas, NumericVector wsigmas){

  //calculating euclidean distance between each observation and each center
  // and idem for the spatial matrix
  NumericMatrix centerdistances(data.nrow(), centers.nrow());
  NumericMatrix Wcenterdistances(data.nrow(), centers.nrow());
  int i;
  for(i = 0 ; i < centers.nrow() ; i++){
    NumericVector ce = centers(i,_);
    centerdistances(_,i) = calcEuclideanDistance2(data, ce) / sigmas(i);
    Wcenterdistances(_,i) = calcEuclideanDistance2(wdata, ce) / wsigmas(i);
  }

  double p = (-1.0/(m-1.0));
  NumericMatrix numerator = power_mat(add_matrices_bycol(centerdistances, (alpha * Wcenterdistances)), p);
  NumericVector denom = rowSums(numerator);

  // finally the total matrix
  NumericMatrix belongmat(numerator.rows(),numerator.cols());
  for(i = 0; i < numerator.cols(); i++){
    NumericVector vals = numerator(_,i) / denom;
    vals[is_na(vals)] = 1;
    belongmat(_,i) = vals;
  }

  return belongmat;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the non spatial version of the classical FCM but with a noisy cluster
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the membership matrix with a noise cluster
//' @name calcBelongMatrixNoisy
//' @description Calculate the membership matrix according to a set of centroids, the observed
//'  data and the fuzziness degree
//' @param centers A matrix or a dataframe representing the centers of the
//'    clusters with p columns and k rows
//' @param data A dataframe or matrix representing the observed data with n rows
//'    and p columns
//' @param m A float representing the fuzziness degree
//' @param delta A float, the value set for delta by the user
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @return A n * k matrix representing the probability of belonging of each
//'    observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcBelongMatrixNoisy(NumericMatrix centers, NumericMatrix data, double m, double delta, NumericVector sigmas){

  int nc = centers.nrow();
  int nd = data.nrow();
  int i,k,j;

  NumericVector delta2(nd, delta);

  // on se constitue une matrice des distances euclidiennes entre chaque
  // observation et chaque groupe
  NumericMatrix distances_to_gp(data.nrow(), nc+1);
  for(i = 0; i < nc ; i++){
    distances_to_gp(_,i) = calcEuclideanDistance2(data, centers(i,_)) / sigmas(i);
  }
  distances_to_gp(_,nc) = delta2;

  NumericMatrix numerator(distances_to_gp.nrow(),distances_to_gp.ncol());
  double p = (1.0/(m-1.0));

  for(i = 0; i < numerator.ncol() ; i++){
    NumericVector ref = distances_to_gp(_,i);
    NumericMatrix base_mat = Rcpp::clone(distances_to_gp);
    for(j = 0; j < base_mat.ncol() ; j++){
      base_mat(_,j) = pow(ref / base_mat(_,j),p);
    }
    numerator(_,i) = rowSums(base_mat);
  }

  NumericMatrix U = power_mat(numerator, -1.0);

  for(i = 0; i < U.ncol(); i ++){
    NumericVector vals = U(_,i);
    vals[is_na(vals)] = 1;
    U(_,i) = vals;
  }

  NumericMatrix U2 = U(_,Range(0,nc-1));
  return(U2);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// For the spatial version of the classical FCM but with a noisy cluster
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title Calculate the membership matrix (spatial version) with a noise cluster
//' @name calcSFCMBelongMatrixNoisy
//' @description Calculate the membership matrix (spatial version) according to a set of
//' centroids, the observed data, the fuzziness degree a neighbouring matrix and
//' a spatial weighting term
//' @param centers A matrix or a dataframe representing the centers of the
//'   clusters with p columns and k rows
//' @param data A matrix representing the observed data with n rows and p columns
//' @param wdata A matrix representing the lagged observed data with n rows and p columns
//' @param m A float representing the fuzziness degree
//' @param alpha A float representing the weight of the space in the analysis (0
//'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
//'   dimensions, 2 is twice the weight for space)
//' @param delta A float, the value set for delta by the user
//' @param sigmas A numeric vector for calculating the robust version of the FCM. Filled with ones
//'    if the classical version is required
//' @param wsigmas Same as sigmas, but calculated on the spatially lagged dataset
//' @return A n * k matrix representing the belonging probabilities of each
//'   observation to each cluster
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix calcSFCMBelongMatrixNoisy(NumericMatrix centers, NumericMatrix data, NumericMatrix wdata, double m, double alpha, double delta, NumericVector sigmas, NumericVector wsigmas){

  //calculating euclidean distance between each observation and each center
  // and idem for the spatial matrix
  int nc = centers.nrow();
  NumericMatrix centerdistances(data.nrow(), nc +1);
  NumericMatrix Wcenterdistances(data.nrow(), nc +1);
  int i;
  for(i = 0 ; i < nc ; i++){
    NumericVector ce = centers(i,_);
    centerdistances(_,i) = calcEuclideanDistance2(data, ce) / sigmas(i);
    Wcenterdistances(_,i) = calcEuclideanDistance2(wdata, ce) / wsigmas(i);
  }

  NumericVector delta2(data.nrow(), delta);
  centerdistances(_,nc) = delta2;
  Wcenterdistances(_,nc) = delta2;

  double p = (-1.0/(m-1.0));
  NumericMatrix numerator = power_mat(add_matrices_bycol(centerdistances, (alpha * Wcenterdistances)), p);
  NumericVector denom = rowSums(numerator);

  // finally the total matrix
  NumericMatrix belongmat(numerator.rows(),numerator.cols());
  for(i = 0; i < numerator.cols(); i++){
    NumericVector vals = numerator(_,i) / denom;
    vals[is_na(vals)] = 1;
    belongmat(_,i) = vals;
  }

  NumericMatrix U2 = belongmat(_,Range(0,nc-1));
  return U2;
}

