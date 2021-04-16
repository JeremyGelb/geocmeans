#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### GFCM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Calculate the generalized membership matrix
#'
#' @description Calculate the generalized membership matrix according to a set of
#' centroids, the observed data, the fuzziness degree, and a beta parameter
#'
#' @param centers A matrix or a dataframe representing the centers of the
#'   clusters with p columns and k rows
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param m A float representing the fuzziness degree
#' @param beta A float for the beta parameter (control speed convergence and classification crispness)
#' @return A n * k matrix representing the belonging probabilities of each
#'   observation to each cluster
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcFGCMBelongMatrix <- function(centers, data, m, beta ){
  #calculating euclidean distance between each observation and each center
  centerdistances <- apply(centers, 1, function(x) {
    return(calcEuclideanDistance(data, x))
  })

  # calculating aj
  aj <- beta * apply(centerdistances,1,min)

  K <- nrow(centers)
  uij <- sapply(1:K, function(i){
    dij <- centerdistances[,i]
    denomVals <- sapply(1:K,function(k){
      dkj <- centerdistances[,k]
      return(((dij-aj) / (dkj-aj))**(1/(m-1)))
    })
    denom <- rowSums(denomVals)
    return(1/denom)
  })
  uij <- ifelse(is.na(uij),1,uij)
  return(uij)
}


#' @title Generalized C-means
#'
#' @description The generalized c-mean algorithm
#'
#' @param data A dataframe with only numerical variable
#' @param k An integer describing the number of cluster to find
#' @param m A float for the fuzziness degree
#' @param beta A float for the beta parameter (control speed convergence and classification crispness)
#' @param maxiter A float for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param standardize A boolean to specify if the variables must be centered and
#'   reduced (default = True)
#' @param verbose A boolean to specify if the messages should be displayed
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @return A named list with :
#'  \itemize{
#'         \item Centers: a dataframe describing the final centers of the groups
#'         \item Belongings: the final membership matrix
#'         \item Groups: a vector with the names of the most likely group for each observation
#'         \item Data: the dataset used to perform the clustering (might be standardized)
#' }
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' result <- GCMeans(dataset,k = 5, m = 1.5, beta = 0.5, standardize = TRUE)
GCMeans <- function(data, k, m, beta, maxiter = 500, tol = 0.01, standardize = TRUE, verbose = TRUE, seed = NULL) {
  # standardize data if required
  if (standardize) {
    if (verbose){
      print("Standardizing the data (set parameter to FALSE to avoid this step)")
    }
    for (i in 1:ncol(data)) {
      data[, i] <- scale(data[, i])
    }
  }

  data <- as.matrix(data)
  results <- main_worker("GFCM", data = data,
                         k = k, m = m, beta = beta,
                         maxiter = maxiter, tol = tol,
                         standardize = standardize, verbose = verbose,
                         seed = seed)
  return(results)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### SGFCM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Calculate the generalized membership matrix (spatial version)
#'
#' @description Calculate the generalized membership matrix (spatial version) according to a set of
#' centroids, the observed data, the fuzziness degree a neighbouring matrix,
#' a spatial weighting term and a beta parameter
#'
#' @param centers A matrix or a dataframe representing the centers of the
#'   clusters with p columns and k rows
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param wdata A dataframe or matrix representing the lagged observed data with
#'   nrows and p columns
#' @param m A float representing the fuzziness degree
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
#' @param beta A float for the beta parameter (control speed convergence and classification crispness)
#' @return A n * k matrix representing the belonging probabilities of each
#'   observation to each cluster
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcSFGCMBelongMatrix <- function(centers, data, wdata, m, alpha, beta ){
  #calculating euclidean distance between each observation and each center
  centerdistances <- apply(centers, 1, function(x) {
    return(calcEuclideanDistance(data, x))
  })

  #same for the lagged data
  Wcenterdistances <- apply(centers, 1, function(x) {
    return(calcEuclideanDistance(wdata, x))
  })

  # calculating aj
  aj <- beta * apply(centerdistances,1,min)

  K <- nrow(centers)
  uij <- sapply(1:K, function(i){
    dij <- centerdistances[,i]
    Wdij <- Wcenterdistances[,i]
    denomVals <- sapply(1:K,function(k){
      dkj <- centerdistances[,k]
      Wdkj <- Wcenterdistances[,k]
      num <- (dij-aj) + alpha * (Wdij)
      denum <- (dkj-aj) + alpha * (Wdkj)
      return(( num / denum)**(1/(m-1)))
    })
    denom <- rowSums(denomVals)
    return(1/denom)
  })
  uij <- ifelse(is.na(uij),1,uij)
  return(uij)
}


#' @title SGFCMeans
#'
#' @description spatial version of the generalized c-mean algorithm (SGFCMeans)
#'
#' @details The implementation is based on the following article : \doi{10.1016/j.dsp.2012.09.016}.\cr
#'
#' the matrix of belonging (u) is calculated as follow \cr
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
#' @param data A dataframe with only numerical variable
#' @param nblistw A list.w object describing the neighbours typically produced
#'   by the spdep package
#' @param k An integer describing the number of cluster to find
#' @param m A float for the fuzziness degree
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
#' @param beta A float for the beta parameter (control speed convergence and classification crispness)
#' @param lag_method A string indicating if a classical lag must be used
#' ("mean") or if a weighted median must be used ("median")
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduced (default = True)
#' @param verbose A boolean to specify if the progress bar should be displayed
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @return A named list with
#' \itemize{
#'         \item Centers: a dataframe describing the final centers of the groups
#'         \item Belongings: the final membership matrix
#'         \item Groups: a vector with the names of the most likely group for each observation
#'         \item Data: the dataset used to perform the clustering (might be standardized)
#' }
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = TRUE)
SGFCMeans <- function(data, nblistw, k, m, alpha, beta, lag_method="mean", maxiter = 500, tol = 0.01, standardize = TRUE, verbose = TRUE, seed = NULL) {

  if(class(nblistw)[[1]] != "listw"){
    stop("the nblistw must be a listw object from spdep package")
  }

  if(lag_method %in% c("mean","median") == FALSE){
    stop("the parameter lag_method must be one of 'mean' or 'median'")
  }


  # standardize data if required
  if (standardize) {
    if (verbose){
      print("Standardizing the data (set parameter to FALSE to avoid this step)")
    }
    for (i in 1:ncol(data)) {
      data[, i] <- scale(data[, i])
    }
  }

  # calculating the lagged dataset
  wdata <- calcLaggedData(data,nblistw,lag_method)

  data <- as.matrix(data)
  wdata <- as.matrix(wdata)

  results <- main_worker('SGFCM', data = data, wdata = wdata,
                         nblistw = nblistw, k = k, m = m, alpha = alpha,
                         beta = beta, lag_method=lag_method, maxiter = maxiter,
                         tol = tol, standardize = standardize, verbose = verbose,
                         seed = seed)
  return(results)
}

