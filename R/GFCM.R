#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### GFCM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Generalized C-means
#'
#' @description The generalized c-mean algorithm
#'
#' @template data_fcm-arg
#' @template k_m-arg
#' @param beta A float for the beta parameter (control speed convergence and classification crispness)
#' @template shared_paramsfcm-arg
#' @template FCMres_return
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' result <- GCMeans(dataset,k = 5, m = 1.5, beta = 0.5, standardize = TRUE)
GCMeans <- function(data, k, m, beta, maxiter = 500, tol = 0.01, standardize = TRUE, verbose = TRUE, init = "random", seed = NULL) {

  #data_class <- class(data)

  if (verbose & standardize){
    print("Standardizing the data (set parameter to FALSE to avoid this step)")
  }

  isRaster <- inherits(data,"list")

  if(isRaster){ # if we have to deal with a raster dataset
    elements <- input_raster_data(dataset  = data,
                                  standardize = standardize)
    data <- elements$data
    missing <- elements$missing
    rst <- elements$rst

  }else if(inherits(data,"data.frame")){ # if we have to deal with a dataframe

    # standardize data if required
    if (standardize) {
      for (i in 1:ncol(data)) {
        data[, i] <- scale(data[, i])
      }
    }
    data <- as.matrix(data)
  }

  results <- main_worker("GFCM", data = data,
                         k = k, m = m, beta = beta,
                         maxiter = maxiter, tol = tol,
                         standardize = standardize, verbose = verbose,
                         seed = seed, init = init)

  # if we are working with rasters, we should set some additional values
  if(isRaster){
    results <- output_raster_data(results, missing, rst)
  }

  return(results)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### SGFCM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title SGFCMeans
#'
#' @description spatial version of the generalized c-mean algorithm (SGFCMeans)
#'
#' @details The implementation is based on the following article : \doi{10.1016/j.dsp.2012.09.016}.\cr
#'
#' the membership matrix (u) is calculated as follow \cr
#' \deqn{u_{ik} = \frac{(||x_{k} - v{_i}||^2 -b_k + \alpha||\bar{x_{k}} - v{_i}||^2)^{(-1/(m-1))}}{\sum_{j=1}^c(||x_{k} - v{_j}||^2 -b_k + \alpha||\bar{x_{k}} - v{_j}||^2)^{(-1/(m-1))}}}
#'
#' the centers of the groups are updated with the following formula
#' \deqn{v_{i} = \frac{\sum_{k=1}^N u_{ik}^m(x_{k} + \alpha\bar{x_{k}})}{(1 + \alpha)\sum_{k=1}^N u_{ik}^m}}
#'
#' with
#' \itemize{
#' \item vi the center of the group vi
#' \item xk the data point k
#' \item xk_bar the spatially lagged data point k
#' \deqn{b_k = \beta \times min(||x_{k} - v||)}
#' }
#'
#' @template data_fcm-arg
#' @template nblistw-arg
#' @template k_m-arg
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
#' @param beta A float for the beta parameter (control speed convergence and classification crispness)
#' @template lag_method-arg
#' @template window-arg
#' @template shared_paramsfcm-arg
#' @template FCMres_return
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = TRUE)
SGFCMeans <- function(data, nblistw = NULL, k, m, alpha, beta, lag_method="mean",
                      window = NULL,
                      maxiter = 500, tol = 0.01,
                      standardize = TRUE, verbose = TRUE, init = "random", seed = NULL) {

  isRaster <- inherits(data, "list")

  if(isRaster==FALSE){
    if(inherits(nblistw, "listw")==FALSE){
      stop("the nblistw must be a listw object from spdep package")
    }
    if(lag_method %in% c("mean","median") == FALSE){
      stop("the parameter lag_method must be one of 'mean' or 'median'")
    }
  }



  if(isRaster){ # if we have to deal with a raster dataset
    elements <- input_raster_data(data,
                                  w = window,
                                  fun = lag_method,
                                  standardize = standardize
    )
    data <- elements$data
    wdata <- elements$wdata
    missing <- elements$missing
    rst <- elements$rst

  }else if(inherits(data,"data.frame")){ # if we have to deal with a dataframe

    # standardize data if required
    if (standardize) {
      for (i in 1:ncol(data)) {
        data[, i] <- scale(data[, i])
      }
    }

    wdata <- calcLaggedData(data,nblistw,lag_method)
    data <- as.matrix(data)
    wdata <- as.matrix(wdata)
  }

  results <- main_worker('SGFCM', data = data, wdata = wdata,
                         nblistw = nblistw, k = k, m = m, alpha = alpha,
                         beta = beta, lag_method=lag_method, maxiter = maxiter,
                         tol = tol, standardize = standardize, verbose = verbose,
                         seed = seed, init = init)

  # if we are working with rasters, we should set some additional values
  if(isRaster){
    results <- output_raster_data(results, missing, rst)
    results$window <- window
  }

  return(results)
}

