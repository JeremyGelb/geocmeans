# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### intermediate functions for cmeans algorithms #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate the euclidean distance between a numeric matrix n * p and a numeric
#' vector of length p
#'
#' @param m A n * p matrix or dataframe with only numeric columns
#' @param v A numeric vector of length p
#' @return A vector of length n giving the euclidean distance between all matrix
#'   row and the vector p
#' @examples
#' #This is an internal function, no example provided
calcEuclideanDistance <- function(m, v) {
    mat1 <- as.matrix(m)
    mat2 <- matrix(as.numeric(v), nrow = nrow(m), ncol = length(v), byrow = TRUE)
    alldistances <- rowSums((mat1 - mat2)^2)
    return(alldistances)
}

#' Intermediar step in the calculation of the belonging matrix
#'
#' @param centerid An integer representing the column to use as reference in the
#'   distance matrix
#' @param alldistances A distance matrix
#' @param m A float, the fuzziness parameter
#' @return m a float representing the fuzzyness degree
#' @examples
#' #This is an internal function, no example provided
calcBelongDenom <- function(centerid, alldistances, m) {
    values <- sapply(1:ncol(alldistances), function(i) {
        div <- (alldistances[, centerid] / alldistances[, i])
        return(div^(2 / (m - 1)))
    })
    return(rowSums(values))
}


#' Calculate the belonging matrix according to a set of centroids, the observed
#' data and the fuzzyness degree
#'
#' @param centers A matrix or a dataframe representing the centers of the
#'   clusters with p columns and k rows
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param m An float representing the fuzzyness degree
#' @return A n * k matrix represening the probability of belonging of each
#'   datapoint to each cluster
#' @examples
#' #This is an internal function, no example provided
calcBelongMatrix <- function(centers, data, m) {
    centerdistances <- as.data.frame(apply(centers, 1, function(x) {
        return(calcEuclideanDistance(data, x))
    }))
    belongingmatrix <- sapply(1:nrow(centers), function(x) {
        return(1 / calcBelongDenom(x, centerdistances, m))
    })
    belongingmatrix[is.na(belongingmatrix)] <- 1
    belongingmatrix[belongingmatrix > 1] <- 1
    return(belongingmatrix)
}


#' Calculate the belonging matrix (spatial version) according to a set of
#' centroids, the observed data, the fuzzyness degree a neighbouring matrix and
#' a spatial ponderation term
#'
#' @param centers A matrix or a dataframe representing the centers of the
#'   clusters with p columns and k rows
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param wdata A dataframe or matrix representing the lagged observed data with
#'   nrows and p columns
#' @param neighbours A list.w object describing the neighbours typically
#'   produced by the spdep package
#' @param m A float representing the fuzzyness degree
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
#' @return A n * k matrix represening the belonging probabilities of each
#'   datapoint to each cluster
#' @examples
#' #This is an internal function, no example provided
calcSFCMBelongMatrix <- function(centers, data, wdata, neighbours, m, alpha) {
    # calculate the original distances (xk-vi)**2
    originaldistances <- apply(centers, 1, function(x) {
        return(calcEuclideanDistance(data, x))
    })
    # calculate the lagged distances (xk_hat-vi)**2
    nbdistances <- apply(centers, 1, function(x) {
        return(calcEuclideanDistance(wdata, x))
    })
    # calculate the denominator
    denom <- rowSums((originaldistances + alpha * nbdistances)^(-1 / (m - 1)))

    # calculate the belonging values
    belongingmatrix <- ((originaldistances + alpha * nbdistances)^(-1 / (m - 1))) / denom

    # check for NA or Inf (if center is located on data point)
    belongingmatrix[is.na(belongingmatrix)] <- 1
    belongingmatrix[belongingmatrix > 1] <- 1
    return(belongingmatrix)
}




#' Calculate the new centroids of the clusters based on the belonging matrix
#'
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param belongmatrix A n X k matrix giving for each observation n, its
#'   probability to belong to the cluster k
#' @param m An integer representing the fuzzyness degree
#' @return A n X k matrix represening the belonging probabilities of each
#'   datapoint to each cluster
#' @examples
#' #This is an internal function, no example provided
calcCentroids <- function(data, belongmatrix, m) {
    centers <- sapply(1:ncol(belongmatrix), function(i) {
        apply(data, 2, function(x) {
            weighted.mean(x, belongmatrix[, i])
        })
    })
    centers <- as.data.frame(t(centers))
    return(centers)
}



#' Calculate the new centroids of the clusters based on the belonging matrix (spatial version)
#'
#' @param data A dataframe or matrix representing the observed data with n rows
#'   and p columns
#' @param wdata A dataframe or matrix representing the lagged observed data with
#'   nrows and p columns
#' @param belongmatrix A n X k matrix giving for each observation n, its
#'   probability to belong to the cluster k
#' @param neighbours A list.w object describing the neighbours typically
#'   produced by the spdep package
#' @param m An integer representing the fuzzyness degree
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
#' @return A n X k matrix represening the belonging probabiblities of each
#'   datapoint to each cluster
#' @examples
#' #This is an internal function, no example provided
calcSWFCCentroids <- function(data, wdata, belongmatrix, neighbours, m, alpha) {
    centers <- sapply(1:ncol(belongmatrix), function(i) {
        weights <- belongmatrix[, i]^m
        sapply(1:ncol(data), function(s) {
            x <- data[, s]
            wx <- wdata[, s]
            v <- (x + alpha * wx)
            return(sum(v * weights) / ((1 + alpha) * sum(weights)))
        })
    })
    centers <- as.data.frame(t(centers))
    return(centers)
}




#' evaluate if the algorithm converged by comparing two successive belongings
#' matrices. Calculate the absolute difference between the matrices and then
#' calculate the mean of each row. If all the values of the final vector are
#' below the fixed tolerance, then return True, else return False
#'
#' @param mat1 A n X k matrix giving for each observation n, its probability to
#'   belong to the cluster k at iteration i
#' @param mat2 A n X k matrix giving for each observation n, its probability to
#'   belong to the cluster k at iteration i+1
#' @param tol a float representing the algorithm tolerance
#' @examples
#' #This is an internal function, no example provided
evaluateMatrices <- function(mat1, mat2, tol) {
    mat1 < -as.matrix(mat1)
    mat2 <- as.matrix(mat2)
    differ <- abs(mat1 - mat2)
    diffobs <- rowSums(differ) / ncol(mat1)
    if (length(diffobs[diffobs >= tol]) > 0) {
        return(FALSE)
    } else (return(TRUE))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Cmeans algorithm function #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' calssical c-mean algorithm
#'
#' @param data A dataframe with only numerical variable
#' @param k An integer describing the number of cluster to find
#' @param m A float for the fuzzyness degree
#' @param maxiter A float for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param standardize A boolean to specify if the variables must be centered and
#'   reduce (default = True)
#' @param verbose A boolean to specify if the messages should be displayed
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @return a named list with :
#'  \itemize{
#'         \item Centers: a dataframe describing the final centers of the groups
#'         \item Belongings: the final belonging matrix
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
    # selecting the original centers from datapoints
    if(is.null(seed)==F){
        set.seed(seed)
    }
    centers <- data[sample(nrow(data), k), ]

    # calculating the first belonging matrix
    belongmatrix <- calcBelongMatrix(centers, data, m)
    CriterioReached <- FALSE

    # starting the loop
    if(verbose){
        pb <- txtProgressBar(1, maxiter, style = 3)
    }
    pb <- txtProgressBar(1, maxiter, style = 3)
    for (i in 1:maxiter) {
        if (verbose){
            setTxtProgressBar(pb, i)
        }
        newcenters <- calcCentroids(data, belongmatrix, m)
        newbelongmatrix <- calcBelongMatrix(newcenters, data, m)
        if (evaluateMatrices(belongmatrix, newbelongmatrix, tol) == FALSE) {
            # if we don't reach convergence criterion
            centers <- newcenters
            belongmatrix <- newbelongmatrix
        } else {
            # if we reach convergence criterion
            if (verbose){
                print("criterion reached")
            }
            CriterioReached <- TRUE
            centers <- newcenters
            belongmatrix <- newbelongmatrix
            break
        }
    }

    if(CriterioReached==FALSE){
        warning("The convergence criterion was not reached within the specified number of steps")
    }
    # calculating the most likely group of earch data point
    DF <- as.data.frame(newbelongmatrix)
    groups <- colnames(DF)[max.col(DF, ties.method = "first")]
    return(list(Centers = centers, Belongings = newbelongmatrix,
                Groups = groups, Data = data))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Spatial Cmeans algorithm function #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' spatial version of the c-mean algorithm (SFCMeans, FCM_S1)
#'
#' The implementation is based on the following article : \url{https://doi.org/10.1016/j.patcog.2006.07.011}.\cr
#' See details for the formulas.
#'
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
#' @param m A float for the fuzzyness degree
#' @param alpha A float representing the weight of the space in the analysis (0
#'   is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
#'   dimensions, 2 is twice the weight for space)
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduce (default = True)
#' @param verbose A boolean to specify if the prossess bar should be displayed
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @return a named list with
#' \itemize{
#'         \item Centers: a dataframe describing the final centers of the groups
#'         \item Belongings: the final belonging matrix
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
SFCMeans <- function(data, nblistw, k, m, alpha, maxiter = 500, tol = 0.01, standardize = TRUE, verbose = TRUE, seed = NULL) {
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
    wdata <- data
    for (Name in names(data)) {
        wdata[[Name]] <- spdep::lag.listw(nblistw, data[[Name]])
    }

    # selecting the original centers from datapoints
    if (is.null(seed)==F){
        set.seed(seed)
    }
    centers <- data[sample(nrow(data), k), ]

    # calculating the first belonging matrix
    belongmatrix <- calcSFCMBelongMatrix(centers, data, wdata, nblistw, m, alpha = alpha)
    CriterioReached <- FALSE

    # starting the loop
    if(verbose){
        pb <- txtProgressBar(1, maxiter, style = 3)
    }
    for (i in 1:maxiter) {
        if(verbose){
            setTxtProgressBar(pb, i)
        }
        newcenters <- calcSWFCCentroids(data, wdata, belongmatrix, nblistw, m, alpha)
        newbelongmatrix <- calcSFCMBelongMatrix(newcenters, data, wdata, nblistw, m, alpha = alpha)
        if (evaluateMatrices(belongmatrix, newbelongmatrix, tol) == FALSE) {
            # if we don't reach convergence criterion
            centers <- newcenters
            belongmatrix <- newbelongmatrix
        } else {
            # if we reach convergence criterion
            if (verbose){
                print("criterion reached")
            }
            CriterioReached <- TRUE
            centers <- newcenters
            break
        }
    }
    # calculating the most likely group of earch data point
    if(CriterioReached==FALSE){
        warning("The convergence criterion was not reached within the specified number of steps")
    }
    DF <- as.data.frame(newbelongmatrix)
    Groups <- colnames(DF)[max.col(DF, ties.method = "first")]
    return(list(Centers = centers, Belongings = newbelongmatrix,
                Groups = Groups, Data = data))
}

