#include "shared_functions.h"
#include <cmath>
using namespace Rcpp;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%% MORAN RELATED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// @title Moran I calculated on a matrix
// @name moranI_matrix
// @param mat a matrix
// @param w the size of the neighbouring window
// @return a double, the value of Moran I
// @keywords internal
// @export
//
// double moranI_matrix(NumericMatrix mat, int w){
//   double xbar = mean(na_omit(mat));
//   //float weight = 1/pow(w*2+1,2);
//   float weight = 1.0;
//   int nr = mat.nrow();
//   int nc = mat.ncol();
//   int r,c,lr,hr,lc,hc;
//   double num = 0.0;
//   double denom = 0.0;
//   double denom2 = 0.0;
//   int N = 0;
//   int wr,wc;
//   NumericMatrix submat, prod;
//   for (r = 0; r < mat.nrow() ; r++){
//     wr = w;
//     lr = r - w;
//     hr = r + w;
//     if(lr < 0){
//       wr = w + lr;
//       lr = 0;
//     }
//     if(hr >= nr){
//       hr = nr-1;
//     }
//     for(c = 0; c < mat.ncol(); c++){
//       wc = w;
//       lc = c - w;
//       hc = c + w;
//       if(lc < 0){
//         wc = w + lc;
//         lc = 0;
//       }
//       if(hc >= nc){
//         hc = nc-1;
//       }
//       float xi =  mat(r,c);
//       if(!Rcpp::traits::is_nan<REALSXP>(xi)){
//         N++;
//         denom = denom + pow((xi - xbar),2);
//         submat = mat(Range(lr,hr),Range(lc,hc));
//         submat(wr,wc) = NA_REAL;
//         prod = (submat - xbar) * (xi - xbar) * weight;
//         num = num + sum(na_omit(prod));
//         denom2 = denom2 + ((submat.ncol() * submat.nrow() -1)* weight);
//
//         }
//       }
//     }
//   double MoranI = (N / denom2) * (num/denom);
//   return MoranI;
// }


//' @title Moran I calculated on a matrix with a given window
//' @name moranI_matrix_window
//' @param mat a matrix
//' @param window the window to use to define neighbours. 0 can be used to indicate that a cell is not a neighbour
//' @return a double, the value of Moran I
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
double moranI_matrix_window(NumericMatrix mat, NumericMatrix window){

  double xbar = mean(na_omit(mat));
  float weight = 1.0;
  int nr = mat.nrow();
  int nc = mat.ncol();
  int r,c,lr,hr,lc,hc;
  double num = 0.0;
  double denom = 0.0;
  double denom2 = 0.0;
  int N = 0;
  int wr,wc;
  int w1 = floor(float(window.nrow())/2.0);
  int w2 = floor(float(window.ncol())/2.0);
  int window_width = window.nrow();
  int window_height = window.ncol();
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;

  NumericMatrix submat, prod, subwindow;
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
      float xi =  mat(r,c);
      if(!Rcpp::traits::is_nan<REALSXP>(xi)){
        N++;
        denom = denom + pow((xi - xbar),2);
        submat = mat(Range(lr,hr),Range(lc,hc));

        subwindow = window(Range(start_window_row,end_window_row),
                           Range(start_window_col,end_window_col));

        //submat(wr,wc) = NA_REAL;
        subwindow(wr,wc) = 0;
        prod = prod_matrices_bycol((submat - xbar) * (xi - xbar),subwindow);
        num = num + sum(na_omit(prod));
        denom2 = denom2 + (sum(subwindow));

      }
    }
  }
  double MoranI = (N / denom2) * (num/denom);
  return MoranI;
}


//' @title Local Moran I calculated on a matrix with a given window
//' @name local_moranI_matrix_window
//' @param mat a matrix
//' @param window the window to use to define neighbours. 0 can be used to indicate that a cell is not a neighbour
//' @return a double, the value of Moran I
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericVector local_moranI_matrix_window(NumericMatrix mat, NumericMatrix window){

  double xbar = mean(na_omit(mat));
  float weight = 1.0;
  int nr = mat.nrow();
  int nc = mat.ncol();
  int r,c,lr,hr,lc,hc;
  int wr,wc;
  int w1 = floor(float(window.nrow())/2.0);
  int w2 = floor(float(window.ncol())/2.0);
  int window_width = window.nrow();
  int window_height = window.ncol();
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;
  NumericVector localI(mat.nrow()*mat.ncol());
  NumericMatrix submat, prod, subwindow;
  int count = 0;

  int n = sum(!is_na(mat));
  //int n = mat.nrow()*mat.ncol();
  double s1 = sum(na_omit(pow(mat-xbar,2))) / float(n-1);
  double comp1, comp2;
  // NOTE : s1 is valid (verified by hand)

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
      float xi =  mat(r,c);
      if(!Rcpp::traits::is_nan<REALSXP>(xi)){
        submat = mat(Range(lr,hr),Range(lc,hc));
        subwindow = window(Range(start_window_row,end_window_row),
                           Range(start_window_col,end_window_col));
        subwindow(wr,wc) = 0;
        comp1 = ((xi - xbar) / s1);
        comp2 = sum(na_omit(prod_matrices_bycol((submat - xbar),subwindow/sum(subwindow))));
        localI(count) = comp1*comp2;
      }else{
        localI(count) = NA_REAL;
      }
      count++;
    }
  }

  return localI;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%% ELSA RELATED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// @title Elsa statistic calculated on a matrix
// @name Elsa_categorical_matrix
// @description method described here : https://doi.org/10.1016/j.spasta.2018.10.001
// @param mat an IntegerMatrix, must be filled with integer, -1 indicates NA values, categories must start at 0
// @param w the size of the neighbouring window
// @param dist a distance matrix between the categories
// @return a NumericVector : the local values of ELSA
// @keywords internal
// @export
//
// NumericVector Elsa_categorical_matrix(IntegerMatrix mat, int w, NumericMatrix dist){
//   int m = dist.nrow();
//   int nr = mat.nrow();
//   int nc = mat.ncol();
//   int r,c,lr,hr,lc,hc,i,j,xi,xj,nn;
//   int N = 0;
//   int wr,wc;
//   double d = max_mat(dist);
//   int counter = 0;
//   double Ea,Ec,pk;
//   NumericVector elsa(mat.nrow()*mat.ncol());
//   IntegerMatrix submat;
//   for (r = 0; r < mat.nrow() ; r++){
//     wr = w;
//     lr = r - w;
//     hr = r + w;
//     if(lr < 0){
//       wr = w + lr;
//       lr = 0;
//     }
//     if(hr >= nr){
//       hr = nr-1;
//     }
//     for(c = 0; c < mat.ncol(); c++){
//       wc = w;
//       lc = c - w;
//       hc = c + w;
//       if(lc < 0){
//         wc = w + lc;
//         lc = 0;
//       }
//       if(hc >= nc){
//         hc = nc-1;
//       }
//       xi =  mat(r,c);
//       if(xi >= 0){
//         Ea = 0;
//         N++;
//         submat = mat(Range(lr,hr),Range(lc,hc));
//         nn = (submat.ncol()*submat.nrow() - 1);
//
//         // now I want to calculate Eci
//         Ec = 0;
//         for(i = 0; i < m; i++){
//           pk = float(sum(submat == i)) / float(nn+1.0);
//           if(pk > 0){
//             Ec = Ec + pk * log2(pk);
//           }
//         }
//
//         if(nn > m){
//           Ec = Ec / log2(float(m));
//         }else{
//           Ec = Ec / log2(float(nn));
//         }
//
//         // now I want to calculate Eai
//         submat(wr,wc) = -1;
//         //here I want to iterate over the neighbours to calculate the sum of their differences
//         //based on the values in dist
//         for(i = 0; i < submat.nrow(); i++){
//           for(j = 0; j < submat.ncol(); j++) {
//             xj = submat(i,j);
//             if(xj >=0 ){
//               Ea = Ea + dist(xi,xj);
//             }
//           }
//         }
//         Ea = Ea / (d * float(nn));
//         elsa(counter) = (-1.0*Ec) * Ea;
//
//       }else{
//       // IF XI is NA
//       elsa(counter) = -1;
//       }
//       counter++;
//     }
//   }
//   return elsa;
//
// }



//' @title Elsa statistic calculated on a matrix with a given window
//' @name Elsa_categorical_matrix_window
//' @description method described here : https://doi.org/10.1016/j.spasta.2018.10.001
//' @param mat an IntegerMatrix, must be filled with integer, -1 indicates NA values, categories must start at 0
//' @param window the window to use to define neighbours. 0 can be used to indicate that a cell is not a neighbour
//' @param dist a distance matrix between the categories
//' @return a NumericVector : the local values of ELSA
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericVector Elsa_categorical_matrix_window(IntegerMatrix mat, IntegerMatrix window, NumericMatrix dist){
  //determiner le centre de la window (width and height of window)
  int window_width = window.nrow();
  int window_height = window.ncol();
  int w1 = floor(float(window.nrow())/2.0);
  int w2 = floor(float(window.ncol())/2.0);
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;

  int m = dist.nrow();
  int nr = mat.nrow();
  int nc = mat.ncol();
  int r,c,lr,hr,lc,hc,i,j,xi,xj,nn;
  int N = 0;
  int wr,wc;
  double d = max_mat(dist);
  int counter = 0;
  double Ea,Ec,pk;
  NumericVector elsa(mat.nrow()*mat.ncol());

  IntegerMatrix submat,subwindow;
  // LogicalMatrix testwindow;
  // iterating over each row
  for (r = 0; r < mat.nrow() ; r++){
    // I need to determine if and how much the window must be cut
    wr = w1;
    start_window_row = 0;
    end_window_row = window_width-1;
    lr = r - w1;
    hr = r + w1;
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
      if(xi >= 0){
        Ea = 0;
        N++;
        submat = mat(Range(lr,hr),Range(lc,hc));
        subwindow = window(Range(start_window_row,end_window_row),
                           Range(start_window_col,end_window_col));
        nn = (submat.ncol()*submat.nrow() - 1);

        //based on that subwindiw, I must update the values in submat
        for(i = 0; i < submat.nrow(); i++){
          for(j = 0; j < submat.ncol(); j++){
              if(subwindow(i,j) == 0){
                submat(i,j) = -1;
              }
          }
        }

        // now I want to calculate Eci
        Ec = 0;
        for(i = 0; i < m; i++){
          pk = float(sum(submat == i)) / float(nn+1.0);
          if(pk > 0){
            Ec = Ec + pk * log2(pk);
          }
        }

        if(nn > m){
          Ec = Ec / log2(float(m));
        }else{
          Ec = Ec / log2(float(nn));
        }

        // now I want to calculate Eai
        submat(wr,wc) = -1;
        //here I want to iterate over the neighbours to calculate the sum of their differences
        //based on the values in dist
        for(i = 0; i < submat.nrow(); i++){
          for(j = 0; j < submat.ncol(); j++) {
            xj = submat(i,j);
            if(xj >=0 ){
              Ea = Ea + dist(xi,xj);
            }
          }
        }
        Ea = Ea / (d * float(nn));
        elsa(counter) = (-1.0*Ec) * Ea;

      }else{
        // IF XI is NA
        elsa(counter) = -1;
      }
      counter++;
    }
  }
  return elsa;
};


//' @title Fuzzy Elsa statistic calculated on a matrix with a given window
//' @name Elsa_fuzzy_matrix_window
//' @description This is an extension to the fuzzy classification case for the Elsa statistic
//' @param mats An array, each slice must contains the membership values of one group
//' @param window the window to use to define neighbours. 0 can be used to indicate that a cell is not a neighbour
//' @param dist a distance matrix between the groups
//' @return a NumericVector : the local values of ELSA
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
NumericVector Elsa_fuzzy_matrix_window(arma::fcube mats, arma::mat window, NumericMatrix dist){
  //determiner le centre de la window (width and height of window)
  int window_width = window.n_rows;
  int window_height = window.n_cols;
  int w1 = floor(float(window_width)/2.0);
  int w2 = floor(float(window_height)/2.0);
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;
  int m = dist.nrow();
  int nr = mats.n_rows;
  int nc = mats.n_cols;
  int r,c,lr,hr,lc,hc,i,j,nn,ii,jj;
  int N = 0;
  int wr,wc;
  double d = max_mat(dist);
  int counter = 0;
  double Ea,Ec,pk,weight,s1;
  NumericVector elsa(nr*nc);

  NumericMatrix submat;
  arma::mat subwindow;
  arma::fvec xi,diftube,xj;
  arma::fcube subcube;
  int L = mats.n_slices;
  // creating an armadillo cube

  // iterating over each row
  for (r = 0; r < nr ; r++){
    // I need to determine if and how much the window must be cut
    wr = w1;
    start_window_row = 0;
    end_window_row = window_width-1;
    lr = r - w1;
    hr = r + w1;
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
      xi =  mats.tube(r,c);

      if(!xi.has_nan()){
        Ea = 0;
        subwindow = window.submat( start_window_row, start_window_col, end_window_row, end_window_col );
        // iterating on each cell in the selected sector
        // and extracting the tube for each cell
        s1 = 0;
        nn = 0;
        subcube = mats.subcube(lr, lc, 0, hr, hc, (L-1));
        for (i = 0; i < subcube.n_rows; i++){
          for(j = 0; j < subcube.n_cols; j++){
            xj = subcube.tube(i,j);
            diftube = abs(xi - xj);
            weight = subwindow(i,j);
            if(weight > 0){
              //calculating the product of each pair and multiply by the dist
              for(ii = 0; ii < xj.n_elem ; ii++){
                for(jj = 0; jj < xj.n_elem ; jj++){
                  s1 += weight * (diftube[ii] * diftube[jj] * dist(ii,jj));
                }
              }
              nn++;
            }
          }
        }

        N++;
        Ea = (s1/2.0) / (float(d)*float(nn-1));

        // now I want to calculate Eci
        Ec = 0;
        for(i = 0; i < m; i++){
          // % stand for element-wise product
          pk = arma::accu((subcube.slice(i) % subwindow)) / float(nn);
          if(pk > 0){
            Ec = Ec + pk * log2(pk);
          }
        }

        if(nn > m){
          Ec = Ec / log2(float(m));
        }else{
          Ec = Ec / log2(float(nn));
        }

        elsa(counter) = (-1.0*Ec) * Ea;

      }else{
        elsa(counter) = NA_REAL;
      }
      counter++;
    }
  }
  return elsa;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%% SPATIAL CONSISTENCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title Adjusted spatial inconsistency index for rasters
//' @name adj_spconsist_arr_window_globstd
//' @param data an arma cube of dimension nr,nc,ns
//' @param memberships an arma cube of dimension nr, nc, ks
//' @param a matrix representing the neighbouring of each pixel
//' @return a double, the adjusted spatial inconsitency index
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
double adj_spconsist_arr_window_globstd(arma::fcube data, arma::fcube memberships, arma::mat window){
  //determiner le centre de la window (width and height of window)
  int window_width = window.n_rows;
  int window_height = window.n_cols;
  int w1 = floor(float(window_width)/2.0);
  int w2 = floor(float(window_height)/2.0);
  int start_window_row = 0, end_window_row = window_width-1, start_window_col = 0, end_window_col = window_height-1;
  int nr = memberships.n_rows;
  int nc = memberships.n_cols;
  int r,c,lr,hr,lc,hc,i,j,nn,ii,jj;
  int N = 0;
  int wr,wc;
  int counter = 0;
  int weight = 0;

  NumericMatrix submat;
  arma::mat subwindow;
  arma::fvec xi,xj,mi,mj;
  arma::fcube subcube_mem, subcube_dat;
  int L = data.n_slices;
  int K = memberships.n_slices;
  double euc_dist;
  double spconsist = 0;

  // **** first serie of iterations to get the overall spatial weight ****
  double tot_weight = 0;
  // iterating over each row
  for (r = 0; r < nr ; r++){
    // I need to determine if and how much the window must be cut
    wr = w1;
    start_window_row = 0;
    end_window_row = window_width-1;
    lr = r - w1;
    hr = r + w1;
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
      xi = data.tube(r,c);

      if(!xi.has_nan()){
        subwindow = window.submat( start_window_row, start_window_col, end_window_row, end_window_col );
        subcube_dat = data.subcube(lr, lc, 0, hr, hc, (L-1));
        // iterating on each cell in the selected sector
        // and extracting the tube for each cell
        for (i = 0; i < subcube_dat.n_rows; i++){
          for(j = 0; j < subcube_dat.n_cols; j++){
            // calculating the euclidean distance between xi and xj
            weight = subwindow(i,j);
            if(weight > 0){
              xj = subcube_dat.tube(i,j);
              if(!xj.has_nan()){
                tot_weight += 1.0 / accu(arma::square(xi - xj));
              }
            }
          }
        }

      }
    }
  }

  // **** second serie of iterations to calculate spconsistency ****


  // iterating over each row
  for (r = 0; r < nr ; r++){
    // I need to determine if and how much the window must be cut
    wr = w1;
    start_window_row = 0;
    end_window_row = window_width-1;
    lr = r - w1;
    hr = r + w1;
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
      xi = data.tube(r,c);
      if(!xi.has_nan()){
        mi =  memberships.tube(r,c);
        subwindow = window.submat( start_window_row, start_window_col, end_window_row, end_window_col );
        subcube_mem = memberships.subcube(lr, lc, 0, hr, hc, (K-1));
        subcube_dat = data.subcube(lr, lc, 0, hr, hc, (L-1));
        int non_zero = accu(subwindow != 0);
        arma::vec dat_dists = zeros(non_zero);
        arma::vec mem_dists = zeros(non_zero);;
        int counter = 0;
        // iterating on each cell in the selected sector
        // and extracting the tube for each cell
        for (i = 0; i < subcube_dat.n_rows; i++){
          for(j = 0; j < subcube_dat.n_cols; j++){
            // calculating the euclidean distance between xi and xj
            weight = subwindow(i,j);
            if(weight > 0){
              xj = subcube_dat.tube(i,j);
              if(!xj.has_nan()){
                mj = subcube_mem.tube(i,j);
                dat_dists(counter) = 1.0 / accu(arma::square(xi - xj));
                mem_dists(counter) = accu(arma::square(mi - mj));
              }else{
                dat_dists(counter) = 0;
                mem_dists(counter) = 0;
              }
              counter++;
            }
          }
        }
        spconsist = spconsist + accu(mem_dists % (dat_dists/tot_weight));
      }
    }
  }
  return spconsist;
}


// mat1 <- rbind(
//     c(0.1,0.5,0.3),
//     c(0.8,0.1,0.1),
//     c(0.7,0.8,0.1)
// )
// mat2 <- 1 - mat1
//
// dat1 <- rbind(
//     c(1,2,3),
//     c(4,5,6),
//     c(1,2,3)
//     )
//
// dat2 <- rbind(
//       c(4,5,6),
//       c(4,4,4),
//       c(2,2,2)
// )
//
// dat3 <- rbind(
//       c(3,3,3),
//       c(1,6,1),
//       c(4,5,7)
//   )
//
// dat_arr <- array(c(dat1,dat2,dat3),c(3,3,3))
// dat_mem <- array(c(mat1,mat2),c(3,3,2))
// window <- matrix(1,nrow = 3, ncol = 3)
// window[2,2] <- 0
// library(raster)
// tiff_location <- system.file("extdata", "Littoral.tif",
//                              package = "geocmeans", mustWork = TRUE)
// full_raster <- brick(tiff_location)
// dataset <- lapply(1:nbands(full_raster), function(i){
//      aband <- raster(tiff_location, band = i)
//      return(aband)
// })
// names(dataset) <- c("blue", "green", "red", "infrared", "SWIR1", "SWIR2")
// FCM_result <- CMeans(dataset, k = 7, m = 1.7, standardize = TRUE, verbose = F, seed = 456, tol = 0.01)
// mats <- lapply(FCM_result$rasters[1:FCM_result$k], raster::as.matrix)
// arr <- array(do.call(c,mats), c(dim(mats[[1]]), length(mats)))
// mats2 <- lapply(dataset, raster::as.matrix)
// arr_dat <-  array(do.call(c,mats2), c(dim(mats2[[1]]), length(mats2)))
// window <- matrix(1,nrow = 3, ncol = 3)
// window[2,2] <- 0
