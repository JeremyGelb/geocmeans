// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calcCentroids
NumericMatrix calcCentroids(NumericMatrix data, NumericMatrix belongmatrix, double m);
RcppExport SEXP _geocmeans_calcCentroids(SEXP dataSEXP, SEXP belongmatrixSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type belongmatrix(belongmatrixSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(calcCentroids(data, belongmatrix, m));
    return rcpp_result_gen;
END_RCPP
}
// calcBelongMatrix
NumericMatrix calcBelongMatrix(NumericMatrix centers, NumericMatrix data, double m);
RcppExport SEXP _geocmeans_calcBelongMatrix(SEXP centersSEXP, SEXP dataSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(calcBelongMatrix(centers, data, m));
    return rcpp_result_gen;
END_RCPP
}
// calcSWFCCentroids
NumericMatrix calcSWFCCentroids(NumericMatrix data, NumericMatrix wdata, NumericMatrix belongmatrix, double m, double alpha);
RcppExport SEXP _geocmeans_calcSWFCCentroids(SEXP dataSEXP, SEXP wdataSEXP, SEXP belongmatrixSEXP, SEXP mSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wdata(wdataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type belongmatrix(belongmatrixSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(calcSWFCCentroids(data, wdata, belongmatrix, m, alpha));
    return rcpp_result_gen;
END_RCPP
}
// calcSFCMBelongMatrix
NumericMatrix calcSFCMBelongMatrix(NumericMatrix centers, NumericMatrix data, NumericMatrix wdata, double m, double alpha);
RcppExport SEXP _geocmeans_calcSFCMBelongMatrix(SEXP centersSEXP, SEXP dataSEXP, SEXP wdataSEXP, SEXP mSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wdata(wdataSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(calcSFCMBelongMatrix(centers, data, wdata, m, alpha));
    return rcpp_result_gen;
END_RCPP
}
// calcFGCMBelongMatrix
NumericMatrix calcFGCMBelongMatrix(NumericMatrix centers, NumericMatrix data, double m, double beta);
RcppExport SEXP _geocmeans_calcFGCMBelongMatrix(SEXP centersSEXP, SEXP dataSEXP, SEXP mSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(calcFGCMBelongMatrix(centers, data, m, beta));
    return rcpp_result_gen;
END_RCPP
}
// calcSFGCMBelongMatrix
NumericMatrix calcSFGCMBelongMatrix(NumericMatrix centers, NumericMatrix data, NumericMatrix wdata, double m, double alpha, double beta);
RcppExport SEXP _geocmeans_calcSFGCMBelongMatrix(SEXP centersSEXP, SEXP dataSEXP, SEXP wdataSEXP, SEXP mSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type wdata(wdataSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(calcSFGCMBelongMatrix(centers, data, wdata, m, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// focal_euclidean_mat_window
NumericMatrix focal_euclidean_mat_window(NumericMatrix mat, NumericMatrix window);
RcppExport SEXP _geocmeans_focal_euclidean_mat_window(SEXP matSEXP, SEXP windowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type window(windowSEXP);
    rcpp_result_gen = Rcpp::wrap(focal_euclidean_mat_window(mat, window));
    return rcpp_result_gen;
END_RCPP
}
// focal_euclidean_list
NumericMatrix focal_euclidean_list(List matrices, NumericMatrix window);
RcppExport SEXP _geocmeans_focal_euclidean_list(SEXP matricesSEXP, SEXP windowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type matrices(matricesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type window(windowSEXP);
    rcpp_result_gen = Rcpp::wrap(focal_euclidean_list(matrices, window));
    return rcpp_result_gen;
END_RCPP
}
// moranI_matrix
double moranI_matrix(NumericMatrix mat, int w);
RcppExport SEXP _geocmeans_moranI_matrix(SEXP matSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(moranI_matrix(mat, w));
    return rcpp_result_gen;
END_RCPP
}
// moranI_matrix_window
double moranI_matrix_window(NumericMatrix mat, NumericMatrix window);
RcppExport SEXP _geocmeans_moranI_matrix_window(SEXP matSEXP, SEXP windowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type window(windowSEXP);
    rcpp_result_gen = Rcpp::wrap(moranI_matrix_window(mat, window));
    return rcpp_result_gen;
END_RCPP
}
// local_moranI_matrix_window
NumericVector local_moranI_matrix_window(NumericMatrix mat, NumericMatrix window);
RcppExport SEXP _geocmeans_local_moranI_matrix_window(SEXP matSEXP, SEXP windowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type window(windowSEXP);
    rcpp_result_gen = Rcpp::wrap(local_moranI_matrix_window(mat, window));
    return rcpp_result_gen;
END_RCPP
}
// Elsa_categorical_matrix
NumericVector Elsa_categorical_matrix(IntegerMatrix mat, int w, NumericMatrix dist);
RcppExport SEXP _geocmeans_Elsa_categorical_matrix(SEXP matSEXP, SEXP wSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dist(distSEXP);
    rcpp_result_gen = Rcpp::wrap(Elsa_categorical_matrix(mat, w, dist));
    return rcpp_result_gen;
END_RCPP
}
// Elsa_categorical_matrix_window
NumericVector Elsa_categorical_matrix_window(IntegerMatrix mat, IntegerMatrix window, NumericMatrix dist);
RcppExport SEXP _geocmeans_Elsa_categorical_matrix_window(SEXP matSEXP, SEXP windowSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type window(windowSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dist(distSEXP);
    rcpp_result_gen = Rcpp::wrap(Elsa_categorical_matrix_window(mat, window, dist));
    return rcpp_result_gen;
END_RCPP
}
// vecmin
double vecmin(NumericVector x);
RcppExport SEXP _geocmeans_vecmin(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecmin(x));
    return rcpp_result_gen;
END_RCPP
}
// vecmax
double vecmax(NumericVector x);
RcppExport SEXP _geocmeans_vecmax(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecmax(x));
    return rcpp_result_gen;
END_RCPP
}
// power_mat
NumericMatrix power_mat(NumericMatrix x, double p);
RcppExport SEXP _geocmeans_power_mat(SEXP xSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(power_mat(x, p));
    return rcpp_result_gen;
END_RCPP
}
// calcEuclideanDistance2
NumericVector calcEuclideanDistance2(NumericMatrix y, NumericVector x);
RcppExport SEXP _geocmeans_calcEuclideanDistance2(SEXP ySEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(calcEuclideanDistance2(y, x));
    return rcpp_result_gen;
END_RCPP
}
// add_matrices_bycol
NumericMatrix add_matrices_bycol(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _geocmeans_add_matrices_bycol(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(add_matrices_bycol(x, y));
    return rcpp_result_gen;
END_RCPP
}
// sub_matrices_bycol
NumericMatrix sub_matrices_bycol(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _geocmeans_sub_matrices_bycol(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(sub_matrices_bycol(x, y));
    return rcpp_result_gen;
END_RCPP
}
// prod_matrices_bycol
NumericMatrix prod_matrices_bycol(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _geocmeans_prod_matrices_bycol(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(prod_matrices_bycol(x, y));
    return rcpp_result_gen;
END_RCPP
}
// div_matrices_bycol
NumericMatrix div_matrices_bycol(NumericMatrix x, NumericMatrix y);
RcppExport SEXP _geocmeans_div_matrices_bycol(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(div_matrices_bycol(x, y));
    return rcpp_result_gen;
END_RCPP
}
// sqrt_matrix_bycol
NumericMatrix sqrt_matrix_bycol(NumericMatrix x);
RcppExport SEXP _geocmeans_sqrt_matrix_bycol(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sqrt_matrix_bycol(x));
    return rcpp_result_gen;
END_RCPP
}
// pow_matrix_bycol
NumericMatrix pow_matrix_bycol(NumericMatrix x, float p);
RcppExport SEXP _geocmeans_pow_matrix_bycol(SEXP xSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< float >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(pow_matrix_bycol(x, p));
    return rcpp_result_gen;
END_RCPP
}
// rowmins_mat
NumericVector rowmins_mat(NumericMatrix x);
RcppExport SEXP _geocmeans_rowmins_mat(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rowmins_mat(x));
    return rcpp_result_gen;
END_RCPP
}
// max_mat
double max_mat(NumericMatrix x);
RcppExport SEXP _geocmeans_max_mat(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(max_mat(x));
    return rcpp_result_gen;
END_RCPP
}
// test_inferior_mat
LogicalMatrix test_inferior_mat(NumericMatrix mat, double t);
RcppExport SEXP _geocmeans_test_inferior_mat(SEXP matSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(test_inferior_mat(mat, t));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_geocmeans_calcCentroids", (DL_FUNC) &_geocmeans_calcCentroids, 3},
    {"_geocmeans_calcBelongMatrix", (DL_FUNC) &_geocmeans_calcBelongMatrix, 3},
    {"_geocmeans_calcSWFCCentroids", (DL_FUNC) &_geocmeans_calcSWFCCentroids, 5},
    {"_geocmeans_calcSFCMBelongMatrix", (DL_FUNC) &_geocmeans_calcSFCMBelongMatrix, 5},
    {"_geocmeans_calcFGCMBelongMatrix", (DL_FUNC) &_geocmeans_calcFGCMBelongMatrix, 4},
    {"_geocmeans_calcSFGCMBelongMatrix", (DL_FUNC) &_geocmeans_calcSFGCMBelongMatrix, 6},
    {"_geocmeans_focal_euclidean_mat_window", (DL_FUNC) &_geocmeans_focal_euclidean_mat_window, 2},
    {"_geocmeans_focal_euclidean_list", (DL_FUNC) &_geocmeans_focal_euclidean_list, 2},
    {"_geocmeans_moranI_matrix", (DL_FUNC) &_geocmeans_moranI_matrix, 2},
    {"_geocmeans_moranI_matrix_window", (DL_FUNC) &_geocmeans_moranI_matrix_window, 2},
    {"_geocmeans_local_moranI_matrix_window", (DL_FUNC) &_geocmeans_local_moranI_matrix_window, 2},
    {"_geocmeans_Elsa_categorical_matrix", (DL_FUNC) &_geocmeans_Elsa_categorical_matrix, 3},
    {"_geocmeans_Elsa_categorical_matrix_window", (DL_FUNC) &_geocmeans_Elsa_categorical_matrix_window, 3},
    {"_geocmeans_vecmin", (DL_FUNC) &_geocmeans_vecmin, 1},
    {"_geocmeans_vecmax", (DL_FUNC) &_geocmeans_vecmax, 1},
    {"_geocmeans_power_mat", (DL_FUNC) &_geocmeans_power_mat, 2},
    {"_geocmeans_calcEuclideanDistance2", (DL_FUNC) &_geocmeans_calcEuclideanDistance2, 2},
    {"_geocmeans_add_matrices_bycol", (DL_FUNC) &_geocmeans_add_matrices_bycol, 2},
    {"_geocmeans_sub_matrices_bycol", (DL_FUNC) &_geocmeans_sub_matrices_bycol, 2},
    {"_geocmeans_prod_matrices_bycol", (DL_FUNC) &_geocmeans_prod_matrices_bycol, 2},
    {"_geocmeans_div_matrices_bycol", (DL_FUNC) &_geocmeans_div_matrices_bycol, 2},
    {"_geocmeans_sqrt_matrix_bycol", (DL_FUNC) &_geocmeans_sqrt_matrix_bycol, 1},
    {"_geocmeans_pow_matrix_bycol", (DL_FUNC) &_geocmeans_pow_matrix_bycol, 2},
    {"_geocmeans_rowmins_mat", (DL_FUNC) &_geocmeans_rowmins_mat, 1},
    {"_geocmeans_max_mat", (DL_FUNC) &_geocmeans_max_mat, 1},
    {"_geocmeans_test_inferior_mat", (DL_FUNC) &_geocmeans_test_inferior_mat, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_geocmeans(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
