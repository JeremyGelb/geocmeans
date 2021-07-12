#include "shared_functions.h"
using namespace Rcpp;


// OLDER VERSION, keeped for debugging purpose
// NumericMatrix focal_euclidean_mat(NumericMatrix mat, int w){
//   int nr = mat.nrow();
//   int nc = mat.ncol();
//   // creating an empty matrix
//   NumericMatrix m2(nr, nc);
//   NumericMatrix submat;
//   // looping over mat values
//   int i,j,lr,hr,lc,hc;
//   for(i = 0; i < nr; ++i) {
//     lr = i - w;
//     hr = i + w;
//     if(lr < 0){
//       lr = 0;
//     }
//     if(hr >= nr){
//       hr = nr-1;
//     }
//     for(j = 0; j < nc; j++){
//       lc = j - w;
//       hc = j + w;
//       if(lc < 0){
//         lc = 0;
//       }
//       if(hc >= nc){
//         hc = nc-1;
//       }
//       submat = mat(Range(lr,hr),Range(lc,hc));
//       m2(i,j) = sum(pow(submat-mat(i,j),2));
//     }
//   }
//   return m2;
// }


//' @title focal euclidean distance on a matrix with a given window
//' @name focal_euclidean_mat_window
//' @param mat a matrix
//' @param window a numeric matrix (squarred)
//' @return a matrix with the euclidean distance of each cell to
//' its neighbours.
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix focal_euclidean_mat_window(NumericMatrix mat, NumericMatrix window){

  int nr = mat.nrow();
  int nc = mat.ncol();
  int r,c,lr,hr,lc,hc;
  int wr,wc;
  int w1 = floor(float(window.nrow())/2.0);
  int w2 = floor(float(window.ncol())/2.0);
  int window_width = window.nrow();
  int window_height = window.ncol();
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;
  NumericMatrix submat, prod, subwindow;
  double xi;
  // creating an empty matrix
  NumericMatrix m2(nr, nc);

  for (r = 0; r < mat.nrow() ; r++){
    wr = w1;
    lr = r - w1;
    hr = r + w1;
    start_window_row = 0;
    end_window_row = window_width-1;
    if(lr < 0){
      wr = w1 + lr;
      start_window_row = lr * -1;
      lr = 0;
    }
    if(hr >= nr){
      end_window_row = window_width - 2 - (hr - nr); // so the width of the window minus excess
      hr = nr-1;
    }
    for(c = 0; c < mat.ncol(); c++){
      start_window_col = 0;
      end_window_col = window_height-1;
      wc = w2;
      lc = c - w2;
      hc = c + w2;
      if(lc < 0){
        wc = w2 + lc;
        start_window_col = lc *-1;
        lc = 0;
      }
      if(hc >= nc){
        end_window_col = window_height - 2 - (hc - nc); // so the height of the window minus excess
        hc = nc-1;
      }
      xi =  mat(r,c);
      if(!Rcpp::traits::is_nan<REALSXP>(xi)){
        submat = mat(Range(lr,hr),Range(lc,hc));
        subwindow = window(Range(start_window_row,end_window_row),
                           Range(start_window_col,end_window_col));
        subwindow(wr,wc) = 0;
        m2(r,c) = sum(na_omit(prod_matrices_bycol(pow_matrix_bycol((submat - xi),2),subwindow)));
      }else{
        m2(r,c) = NA_REAL;
      }
    }
  }
  return m2;
}


//' @title focal euclidean distance on a list of matrices
//' @name focal_euclidean
//' @param matrices a List of matrices with the same dimensions
//' @param window a numeric matrix
//' @return a matrix with the euclidean distance of each cell to
//' its neighbours.
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix focal_euclidean_list(List matrices, NumericMatrix window){
  NumericMatrix mat = matrices[0];
  int nr = mat.nrow();
  int nc = mat.ncol();
  // creating an empty matrix
  NumericMatrix m2(nr, nc);

  int i;
  for(i = 0; i< matrices.length(); i++){
    NumericMatrix m1 = matrices[i];
    m2 = add_matrices_bycol(m2, focal_euclidean_mat_window(m1,window));
  }
  return m2;
}






