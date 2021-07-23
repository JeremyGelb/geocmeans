#include "geocmeans.h"

//' @title minimum of a vector
//' @name vecmin
//' @param x a NumericVector
//' @return a double
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
double vecmin(NumericVector x){
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}


//' @title maximum of a vector
//' @name vecmin
//' @param x a NumericVector
//' @return a double
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
double vecmax(NumericVector x){
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}



//' @title power of a matrix
//' @name power_mat
//' @param x a matrix
//' @param p a float
//' @return x ** p
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
NumericMatrix power_mat(NumericMatrix x, double p){
  int nr = x.nrow();
  int nc = x.ncol();
  NumericMatrix xp(nr, nc);
  int i;
  for(i = 0; i < nc; i++){
    xp(_,i) = pow(x(_,i),p);
  }
  return xp;
}


//' @title euclidean distance between rows of a matrix and a vector
//' @name calcEuclideanDistance2
//' @param y a matrix
//' @param x a vector (same length as ncol(matrix))
//' @return a vector (same length as nrow(matrix))
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
NumericVector calcEuclideanDistance2(NumericMatrix y, NumericVector x){
  int nr = y.nrow();
  int i;
  NumericVector result(nr);
  for(i = 0; i < nr ; i++){
    result(i) = sum(pow(x - y(i,_), 2));
  }
  return result;
}

//' @title euclidean distance between rows of a matrix and a vector (arma mode)
//' @name calcEuclideanDistance3
//' @param y a matrix
//' @param x a vector (same length as ncol(matrix))
//' @return a vector (same length as nrow(matrix))
//' @export
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat calcEuclideanDistance3(arma::mat y, arma::mat x){
  return arma::sum(arma::pow(y.each_row() - x,2),1);
}



//' @title sum of two matrices by column
//' @name add_matrices_bycol
//' @param x a matrix
//' @param y a matrix with the same dimensions
//' @return a matrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix add_matrices_bycol(NumericMatrix x, NumericMatrix y){
  NumericMatrix z(x.nrow(), y.ncol());
  int i;
  for(i = 0; i < x.ncol() ; i ++){
    z(_,i) = x(_,i) + y(_,i);
  }
  return z;
}

//' @title substraction of two matrices by column
//' @name sub_matrices_bycol
//' @param x a matrix
//' @param y a matrix with the same dimensions
//' @return a matrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix sub_matrices_bycol(NumericMatrix x, NumericMatrix y){
  NumericMatrix z(x.nrow(), y.ncol());
  int i;
  for(i = 0; i < x.ncol() ; i ++){
    z(_,i) = x(_,i) - y(_,i);
  }
  return z;
}



//' @title element wise product of two matrices by column
//' @name prod_matrices_bycol
//' @param x a matrix
//' @param y a matrix with the same dimensions
//' @return a matrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix prod_matrices_bycol(NumericMatrix x, NumericMatrix y){
  NumericMatrix z(x.nrow(), y.ncol());
  int i;
  for(i = 0; i < x.ncol() ; i ++){
    z(_,i) = x(_,i) * y(_,i);
  }
  return z;
}


//' @title element wise division of two matrices by column
//' @name div_matrices_bycol
//' @param x a matrix
//' @param y a matrix with the same dimensions
//' @return a matrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix div_matrices_bycol(NumericMatrix x, NumericMatrix y){
  NumericMatrix z(x.nrow(), y.ncol());
  int i;
  for(i = 0; i < x.ncol() ; i ++){
    z(_,i) = x(_,i) / y(_,i);
  }
  return z;
}


//' @title element wise square root of a matrix by column
//' @name sqrt_matrix_bycol
//' @param x a matrix
//' @return a matrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix sqrt_matrix_bycol(NumericMatrix x){
  NumericMatrix z(x.nrow(), x.ncol());
  int i;
  for(i = 0; i < x.ncol() ; i ++){
    z(_,i) = sqrt(x(_,i));
  }
  return z;
}


//' @title element wise power of a matrix by column
//' @name pow_matrices_bycol
//' @param x a matrix
//' @param p the exponent
//' @return a matrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix pow_matrix_bycol(NumericMatrix x, float p){
  NumericMatrix z(x.nrow(), x.ncol());
  int i;
  for(i = 0; i < x.ncol() ; i ++){
    z(_,i) = pow(x(_,i),p);
  }
  return z;
}





//' @title minimum of each row of a matrix
//' @name rowmins_mat
//' @param x a matrix
//' @return a NumericVector
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericVector rowmins_mat(NumericMatrix x){
  int n = x.nrow();
  NumericVector mins(n);
  int i;
  for (i = 0; i < n; i++){
    mins(i)=(vecmin(x(i,_)));
  }
  return mins;
}



//' @title maximum in a matrix
//' @name max_mat
//' @param x a matrix
//' @return a double
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
double max_mat(NumericMatrix x){
  int n = x.ncol();
  NumericVector col_maxs(n);
  int i;
  for (i = 1; i < n; i++){
    col_maxs(i)=(vecmax(x(_,i)));
  }
  return vecmax(col_maxs);
}


//' @title create a logical matrix with inferior comparison
//' @name test_inferior_mat
//' @param mat a matrix
//' @param t a double to compare
//' @return a LogicalMatrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
LogicalMatrix test_inferior_mat(NumericMatrix mat, double t){
  int n = mat.nrow(), m = mat.ncol();
  LogicalMatrix result(n, m);
  for ( int i = 0; i < n; ++i ) {
    for ( int j = 0; j < m; ++j ) {
      result(i, j) = mat(i, j) < t;
    }
  }
  return result;
}

//' @title create a matrix by multiplying a vector by its elements one by one as rows
//' @name vector_out_prod
//' @param x a vector
//' @return a NumericMatrix
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix vector_out_prod(NumericVector x){
  NumericMatrix res(x.length(),x.length());
  for(int i = 0; i < x.length(); i++){
      res(i,_) = x(i) * x;
  }
  return res;
}
