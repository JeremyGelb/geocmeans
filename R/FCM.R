# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Cmeans algorithm function #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title C-means
#'
#' @description The classical c-mean algorithm
#'
#' @template data_fcm-arg
#' @template k_m-arg
#' @template shared_paramsfcm-arg
#' @template FCMres_return
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' result <- CMeans(dataset,k = 5, m = 1.5, standardize = TRUE)
CMeans <- function(data, k, m, maxiter = 500, tol = 0.01, standardize = TRUE,
                   robust = FALSE,
                   noise_cluster = FALSE, delta = NULL,
                   verbose = TRUE, init = "random", seed = NULL) {

    #data_class <- class(data)

    if (verbose & standardize){
        print("Standardizing the data (set parameter to FALSE to avoid this step)")
    }

    isRaster <- inherits(data, "list")

    if(isRaster){ # if we have to deal with a raster dataset
        elements <- input_raster_data(dataset  = data,
                                      standardize = standardize)
        data <- elements$data
        missing <- elements$missing
        rst <- elements$rst

    }else if(inherits(data, "data.frame")){ # if we have to deal with a dataframe

        # standardize data if required
        if (standardize) {
            for (i in 1:ncol(data)) {
                data[, i] <- scale(data[, i])
            }
        }
        data <- as.matrix(data)
    }

    results <- main_worker("FCM", data = data, k = k, m = m,
                           maxiter = maxiter, tol = tol, standardize = standardize,
                           verbose = verbose, seed = seed, init = init,
                           robust = robust, noise_cluster = noise_cluster,
                           delta = delta
                           )

    # if we are working with rasters, we should set some additional values
    if(isRaster){
        results <- output_raster_data(results, missing, rst)
    }

    return(results)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Spatial Cmeans algorithm function #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title SFCMeans
#'
#' @description spatial version of the c-mean algorithm (SFCMeans, FCM_S1)
#'
#' @details The implementation is based on the following article : \doi{10.1016/j.patcog.2006.07.011}.\cr
#'
#' the membership matrix (u) is calculated as follow \cr
#' \deqn{u_{ik} = \frac{(||x_{k} - v{_i}||^2 + \alpha||\bar{x_{k}} - v{_i}||^2)^{(-1/(m-1))}}{\sum_{j=1}^c(||x_{k} - v{_j}||^2 + \alpha||\bar{x_{k}} - v{_j}||^2)^{(-1/(m-1))}}}
#'
#' the centers of the groups are updated with the following formula
#' \deqn{v_{i} = \frac{\sum_{k=1}^N u_{ik}^m(x_{k} + \alpha\bar{x_{k}})}{(1 + \alpha)\sum_{k=1}^N u_{ik}^m}}
#'
#' with
#' \itemize{
#' \item vi the center of the group vi
#' \item xk the data point k
#' \item xk_bar the spatially lagged data point k
#' }
#'
#' @template data_fcm-arg
#' @template nblistw-arg
#' @template k_m-arg
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
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
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
SFCMeans <- function(data, nblistw = NULL, k, m, alpha, lag_method="mean",
                     window = NULL,
                     noise_cluster = FALSE, delta = NULL,
                     maxiter = 500, tol = 0.01, standardize = TRUE, robust = FALSE,
                     verbose = TRUE, init = "random", seed = NULL) {

    isRaster <- inherits(data, "list")
    if(isRaster == FALSE){
        if(inherits(nblistw, "listw") == FALSE ){
            stop("the nblistw must be a listw object from spdep package")
        }

        if(lag_method %in% c("mean","median") == FALSE){
            stop("the parameter lag_method must be one of 'mean' or 'median'")
        }
    }

    if (verbose & standardize){
        print("Standardizing the data (set parameter to FALSE to avoid this step)")
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

    }else if(inherits(data, "data.frame")){ # if we have to deal with a dataframe

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
    # calculating the final results
    results <- main_worker("SFCM", data = data, wdata = wdata,
                           nblistw = nblistw, k = k, m = m, alpha = alpha,
                           lag_method=lag_method, maxiter = maxiter, tol = tol,
                           standardize = standardize, verbose = verbose,
                           robust = robust,
                           noise_cluster = noise_cluster, delta = delta,
                           seed = seed, init = init)

    # if we are working with rasters, we should set some additional values
    if(isRaster){
        results <- output_raster_data(results, missing, rst)
        results$window <- window
    }

    return(results)
}

