#include "geocmeans.h"

#ifndef shared_functions
#define shared_functions

// #include <Rcpp.h>
// using namespace Rcpp;
// using namespace arma;

NumericMatrix power_mat(NumericMatrix x, double p);

NumericVector calcEuclideanDistance2(NumericMatrix y, NumericVector x);

arma::vec calcEuclideanDistance3(arma::mat y, arma::mat x);

NumericMatrix add_matrices_bycol(NumericMatrix x, NumericMatrix y);

NumericMatrix sub_matrices_bycol(NumericMatrix x, NumericMatrix y);

NumericMatrix prod_matrices_bycol(NumericMatrix x, NumericMatrix y);

NumericMatrix div_matrices_bycol(NumericMatrix x, NumericMatrix y);

NumericMatrix sqrt_matrix_bycol(NumericMatrix x);

NumericMatrix pow_matrix_bycol(NumericMatrix x, float p);

NumericVector rowmins_mat(NumericMatrix x);

LogicalMatrix test_inferior_mat(NumericMatrix mat, double t);

NumericMatrix vector_out_prod(NumericVector x);

double calcDeltaNoiseClust(NumericMatrix data, NumericMatrix centers, double lambda);

double max_mat(NumericMatrix x);

double vecmin(NumericVector x);

double vecmax(NumericVector x);

NumericVector calcRobustSigmas(NumericMatrix data, NumericMatrix belongmatrix, NumericMatrix centers, double m);

#endif
