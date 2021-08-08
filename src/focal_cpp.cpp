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
//' @param window a numeric matrix (squared)
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



//' @title focal euclidean distance on a matrix with a given window for a cube
//' @name focal_euclidean_arr_window
//' @param mat an array (cube)
//' @param window a numeric matrix (squared)
//' @return a matrix with the euclidean distance of each cell to
//' its neighbours.
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericMatrix focal_euclidean_arr_window(arma::cube mat, arma::mat window){

  int nr = mat.n_rows;
  int nc = mat.n_cols;
  int nl = mat.n_slices;
  int r,c,lr,hr,lc,hc,i,j;
  int wr,wc;
  int w1 = floor(float(window.n_rows)/2.0);
  int w2 = floor(float(window.n_cols)/2.0);
  int window_width = window.n_rows;
  int window_height = window.n_cols;
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;
  arma::mat subwindow;
  double distsums;
  arma::vec xi,xj;
  arma::cube subcube;
  // creating an empty matrix
  NumericMatrix m2(nr, nc);

  for (r = 0; r < nr ; r++){
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
    for(c = 0; c < nc; c++){
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
      xi =  mat.tube(r,c);
      if(!xi.has_nan()){
        subcube = mat.subcube(lr, lc, 0, hr, hc, (nl-1));
        distsums = 0;
        subwindow = window.submat( start_window_row, start_window_col, end_window_row, end_window_col);
        subwindow(wr,wc) = 0;

        for(i = 0; i < subcube.n_rows; i++){
          for(j = 0; j < subcube.n_cols; j++){
            xj = subcube.tube(i,j);
            distsums = distsums + accu(pow(xi - xj,2)) * subwindow(i,j);
          }
        }

        m2(r,c) = distsums;
      }else{
        m2(r,c) = NA_REAL;
      }
    }
  }
  return m2;
}


//' @title focal mean weighted by inverse of euclidean distance on a cube
//' @name focal_adj_mean_arr_window
//' @param mat an array (cube)
//' @param window a numeric matrix (squared)
//' @return a lagged version of the original cube
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
arma::cube focal_adj_mean_arr_window(arma::cube mat, arma::mat window){

  int nr = mat.n_rows;
  int nc = mat.n_cols;
  int nl = mat.n_slices;
  int r,c,lr,hr,lc,hc,i,j;
  int wr,wc;
  int w1 = floor(float(window.n_rows)/2.0);
  int w2 = floor(float(window.n_cols)/2.0);
  int window_width = window.n_rows;
  int window_height = window.n_cols;
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;
  arma::mat subwindow;
  double distsums;
  arma::vec xi,xj;
  arma::cube subcube;
  // creating the output cube with the lagged version of the data
  arma::cube out_cube(mat.n_rows,mat.n_cols,mat.n_slices);

  for (r = 0; r < nr ; r++){
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
    for(c = 0; c < nc; c++){
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
      xi =  mat.tube(r,c);
      if(!xi.has_nan()){
        subcube = mat.subcube(lr, lc, 0, hr, hc, (nl-1));
        subwindow = window.submat(start_window_row, start_window_col, end_window_row, end_window_col);
        arma::mat weightwindow(subwindow.n_rows, subwindow.n_cols);
        subwindow(wr,wc) = 0;
        // calculating the weight window
        for(i = 0; i < subcube.n_rows; i++){
          for(j = 0; j < subcube.n_cols; j++){
            xj = subcube.tube(i,j);
            double dist = accu(pow(xi - xj,2)) * subwindow(i,j);
            if(dist == 0){
              // NOTE : this extrem case is not really satisfactory
              weightwindow(i,j) = 1.0 / 0.00000001;
            }else{
              weightwindow(i,j) = 1.0/dist;
            }
          }
        }
        weightwindow(wr,wc) = 0;
        // standardisation of the weights
        weightwindow = weightwindow / accu(weightwindow);
        // and now for each slices, calculate the value (mean)
        for(i = 0; i < subcube.n_slices; i++){
          out_cube(r,c,i) = accu(subcube.slice(i) % weightwindow % subwindow);
        }
      }else{
        for(i = 0; i < mat.n_slices; i++){
          out_cube(r,c,i) = NA_REAL;
        }
      }
    }
  }
  return out_cube;
}



