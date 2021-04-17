# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Cmeans algorithm function #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Calculate the membership matrix
#'
#' @description Calculate the membership matrix according to a set of centroids, the observed
#' data and the fuzziness degree
#'
#' @param centers A matrix or a dataframe representing the centers of the
#'   clusters with p columns and k rows
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param m A float representing the fuzziness degree
#' @return A n * k matrix representing the probability of belonging of each
#'   observation to each cluster
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
#'
calcBelongMatrix <- function(centers, data, m) {
    #calculating euclidean distance between each observation and each center
    centerdistances <- apply(centers, 1, function(x) {
        return(calcEuclideanDistance(data, x))
    })
    numerator <- centerdistances**(-1/(m-1))
    #calculating the denominator
    denom <- rowSums(numerator)
    #finally the total matrix
    belongmat <- numerator / denom
    # check for NA or Inf (if center is located on data point)
    belongmat[is.na(belongmat)] <- 1
    belongmat[belongmat > 1] <- 1
    return(belongmat)
}


#' @title Calculate the centroids
#'
#' @description Calculate the new centroids of the clusters based on the membership matrix
#'
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param belongmatrix A n X k matrix giving for each observation n, its
#'   probability to belong to the cluster k
#' @param m An integer representing the fuzziness degree
#' @return A n X k matrix representing the belonging probabilities of each
#'   observation to each cluster
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcCentroids <- function(data, belongmatrix, m){
    centers <- sapply(1:ncol(belongmatrix),function(i){
        apply(data,2,function(x){
            sum(x *  belongmatrix[, i]**m) / sum(belongmatrix[, i]**m)
        })
    })
    centers <- as.data.frame(t(centers))
    return(centers)
}


#' @title C-means
#'
#' @description The clasical c-mean algorithm
#'
#' @param data A dataframe with only numerical variable
#' @param k An integer describing the number of cluster to find
#' @param m A float for the fuzziness degree
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
#' result <- CMeans(dataset,k = 5, m = 1.5, standardize = TRUE)
CMeans <- function(data, k, m, maxiter = 500, tol = 0.01, standardize = TRUE, verbose = TRUE, seed = NULL) {
    # standardize data if required
    if (standardize) {
        if(verbose){
            print("Standardizing the data (set parameter to FALSE to avoid this step)")
        }
        for (i in 1:ncol(data)) {
            data[, i] <- scale(data[, i])
        }
    }
    data <- as.matrix(data)
    results <- main_worker("FCM", data = data, k = k, m = m,
                           maxiter = maxiter, tol = tol, standardize = standardize,
                           verbose = verbose, seed = seed
                           )
    return(results)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Spatial Cmeans algorithm function #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Calculate the membership matrix (spatial version)
#'
#' @description Calculate the membership matrix (spatial version) according to a set of
#' centroids, the observed data, the fuzziness degree a neighbouring matrix and
#' a spatial weighting term
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
#' @return A n * k matrix representing the belonging probabilities of each
#'   observation to each cluster
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcSFCMBelongMatrix <- function(centers, data, wdata, m, alpha) {
    #calculating euclidean distance between each observation and each center
    centerdistances <- apply(centers, 1, function(x) {
        return(calcEuclideanDistance(data, x))
    })
    # and now for the lagged data
    Wcenterdistances <- apply(centers, 1, function(x) {
        return(calcEuclideanDistance(wdata, x))
    })

    numerator <- (centerdistances + alpha * Wcenterdistances)**(-1/(m-1))
    #calculating the denominator
    denom <- rowSums(numerator)
    #finally the total matrix
    belongmat <- numerator / denom

    # check for NA or Inf (if center is located on data point)
    belongmat[is.na(belongmat)] <- 1
    belongmat[belongmat > 1] <- 1
    return(belongmat)
}


#' @title Calculate the centroids (spatial version)
#'
#' @description Calculate the new centroids of the clusters based on the membership matrix (spatial version)
#'
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param wdata A dataframe or matrix representing the lagged observed data with
#'   nrows and p columns
#' @param belongmatrix A n X k matrix giving for each observation n, its
#'   probability to belong to the cluster k
#' @param m An integer representing the fuzziness degree
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
#' @return A n X k matrix representing the belonging probabilities of each
#'   observation to each cluster
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcSWFCCentroids <- function(data, wdata, belongmatrix, m, alpha) {
    powered <- belongmatrix**m
    wdata_alpha <- alpha*wdata
    centers <- sapply(1:ncol(belongmatrix),function(i){
        sapply(1:ncol(data),function(y){
            x <- data[,y]
            wx <- wdata_alpha[,y]
            s1 <- powered[, i]
            numerator <- sum( (x + wx) *  s1)
            denom <- ( (1+alpha) * sum(s1))
            return(numerator/denom)
        })
    })
    centers <- as.data.frame(t(centers))
    return(centers)
}



#' @title SFCMeans
#'
#' @description spatial version of the c-mean algorithm (SFCMeans, FCM_S1)
#'
#' @details The implementation is based on the following article : \doi{10.1016/j.patcog.2006.07.011}.\cr
#'
#' the matrix of belonging (u) is calculated as follow \cr
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
#' @param data A dataframe with only numerical variable
#' @param nblistw A list.w object describing the neighbours typically produced
#'   by the spdep package
#' @param k An integer describing the number of cluster to find
#' @param m A float for the fuzziness degree
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
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
#'         \item Belongings: the final bmembership matrix
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
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
SFCMeans <- function(data, nblistw, k, m, alpha, lag_method="mean", maxiter = 500, tol = 0.01, standardize = TRUE, verbose = TRUE, seed = NULL) {

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
    results <- main_worker("SFCM", data = data, wdata = wdata,
                           nblistw = nblistw, k = k, m = m, alpha = alpha,
                           lag_method=lag_method, maxiter = maxiter, tol = tol,
                           standardize = standardize, verbose = verbose, seed = seed)
    return(results)
}

