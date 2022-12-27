#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### intermediate general functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Lagged Data
#'
#' @description Calculate Wx, the spatially lagged version of x, by a neighbouring matrix W.
#'
#' @param x A dataframe with only numeric columns
#' @param nblistw The listw object (spdep like) used to calculate WY
#' @param method A string indicating if a classical lag must be used
#' ("mean") or if a weighted median must be used ("median")
#' @return A lagged version of x
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcLaggedData <- function(x,nblistw,method="mean"){
  if (method=="mean"){
    wx <- x
    for (Name in names(x)) {
      wx[[Name]] <- spdep::lag.listw(nblistw, x[[Name]])
    }
    return(wx)
  }else if (method=="median"){
    nb <- nblistw$neighbours
    weights <- nblistw$weights
    all_values <- lapply(1:length(nblistw$neighbours),function(i){
      ids <- nb[[i]]
      w <- weights[[i]]
      obs <- x[ids,]
      if (ncol(x) == 1){
        values <- reldist::wtd.quantile(obs,q=0.5,weight=w)
      }else{
        values <- apply(obs,2,function(column){
          return(reldist::wtd.quantile(column,q=0.5,weight=w))
        })
      }
      return(values)
    })
    wx <- data.frame(do.call(rbind,all_values))
    return(wx)
  }else {
    stop("The method used to calculate lagged values must be mean or median")
  }
}



#' @title Calculate the Euclidean distance
#'
#' @description Calculate the euclidean distance between a numeric matrix n * p and a numeric
#' vector of length p
#' @param m A n * p matrix or dataframe with only numeric columns
#' @param v A numeric vector of length p
#' @return A vector of length n giving the euclidean distance between all matrix
#'   row and the vector p
#' @importFrom matrixStats colSums2
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcEuclideanDistance <- function(m, v) {
  v <- as.numeric(v)
  alldistances <- colSums2((t(m)-v)**2)
  return(alldistances)
}


#' @title Matrix evaluation
#'
#' @description Evaluate if the algorithm converged by comparing two successive membership
#' matrices. Calculate the absolute difference between the matrices and then
#' calculate the max of each row. If all the values of the final vector are
#' below the fixed tolerance, then return True, else return False
#'
#' @param mat1 A n X k matrix giving for each observation n, its probability to
#'   belong to the cluster k at iteration i
#' @param mat2 A n X k matrix giving for each observation n, its probability to
#'   belong to the cluster k at iteration i+1
#' @param tol A float representing the algorithm tolerance
#' @return A boolean, TRUE if the test is passed, FALSE otherwise
#' @keywords internal
#' @importFrom matrixStats rowMaxs
#' @examples
#' #This is an internal function, no example provided
evaluateMatrices <- function(mat1, mat2, tol) { #nocov start
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  differ <- abs(mat1 - mat2)
  diffobs <- rowMaxs(differ)
  if (length(diffobs[diffobs >= tol]) > 0) {
    return(FALSE)
  } else{
    return(TRUE)
  }
} #nocov end



#' @title kpp centers selection
#'
#' @description Select the initial centers of centroids by using the k++ approach
#' as suggested in this article: http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf
#'
#' @param data The dataset used in the classification
#' @param k The number of groups for the classification
#' @return a DataFrame, each row is the center of a cluster
#' @keywords internal
#' @importFrom matrixStats rowMins
#' @examples
#' #This is an internal function, no example provided
kppCenters <- function(data,k){
  # step1: select the first center at random
  centers <- t(as.matrix(data[sample(1:nrow(data),size = 1),], nrow = 1))
  # select the other centers
  for (i in 2:k){
    dists <- apply(centers,1,function(ci){
      return(calcEuclideanDistance2(data,ci))
    })
    Dx <- rowMins(dists)
    probs <- (Dx**2) / sum((Dx**2))
    new_center <- data[sample(1:nrow(data),size = 1,prob = probs),]
    centers <- rbind(centers,new_center)
  }
  return(centers)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Main worker function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Main worker function
#'
#' @description Execution of the classification algorithm
#'
#' @param algo A string indicating the algorithm to use (one of FCM, GFCM, SGFCM)
#' @param ... all the required arguments for the algorithm to use
#' @keywords internal
#' @return A named list with
#' \itemize{
#'         \item Centers: a dataframe describing the final centers of the groups
#'         \item Belongings: the final membership matrix
#'         \item Groups: a vector with the names of the most likely group for each observation
#'         \item Data: the dataset used to perform the clustering (might be standardized)
#' }
#' @examples
#' #This is an internal function, no example provided
main_worker <- function(algo, ...){
  dots <- list(...)
  verbose <- dots$verbose
  data <- dots$data
  dots$data <- NULL
  maxiter <- dots$maxiter
  tol <- dots$tol
  k <- dots$k
  update_belongs <- NULL
  udpdate_centers <- NULL
  robust <- dots$robust
  dots$sigmas <- rep(1, k)
  dots$wsigmas <- rep(1, k)

  if(is.null(dots$noise_cluster)){
    dots$noise_cluster <- FALSE
  }

  #checking if the parameters are ok
  sanity_check(dots,data)

  if(is.null(dots$delta)){
    dots$delta <- -1
  }


  #selecting the functions for calculus
  if(algo == "FCM"){
    update_belongs <- belongsFCM
    udpdate_centers <- centersFCM
    params <- list(
      k = k,
      m = dots$m,
      algo = "FCM"
    )

  }else if(algo=="GFCM"){
    update_belongs <- belongsGFCM
    udpdate_centers <- centersGFCM
    params <- list(
      k = k,
      m = dots$m,
      beta = dots$beta,
      algo = "GFCM"
    )

  }else if (algo=="SFCM"){
    update_belongs <- belongsSFCM
    udpdate_centers <- centersSFCM
    params <- list(
      k = k,
      m = dots$m,
      alpha = dots$alpha,
      nblistw = dots$nblistw,
      lag_method = dots$lag_method,
      algo = "SFCM"
    )

  }else if (algo=="SGFCM"){
    update_belongs <- belongsSGFCM
    udpdate_centers <- centersSGFCM
    params <- list(
      k = k,
      m = dots$m,
      alpha = dots$alpha,
      nblistw = dots$nblistw,
      beta = dots$beta,
      lag_method = dots$lag_method,
      algo = "SGFCM"
    )
  }

  params$robust <- robust

  # selecting the original centers from observations
  if (is.null(dots$seed)==FALSE){
    set.seed(dots$seed)
  }

  if(dots$init == "random"){
    centers <- data[sample(nrow(data), k), ]
  }else if (dots$init == "kpp") {
    centers <- kppCenters(data,k)
  }

  if(is.null(dim(centers))){
    centers <- matrix(centers, ncol = 1)
  }

  # calculating the first membership matrix
  belongmatrix <- update_belongs(data, centers, dots)
  CriterioReached <- FALSE

  # starting the loop
  if(verbose){
    pb <- txtProgressBar(1, maxiter, style = 3)
  }
  for (i in 1:maxiter) {
    if(verbose){
      setTxtProgressBar(pb, i)
    }
    newcenters <- udpdate_centers(data, centers, belongmatrix, dots)

    # for the robust version, we need to calculate here the sigmas
    if(robust){
      dots$sigmas <- calcRobustSigmas(data, belongmatrix, centers, m = dots$m)
      if(is.null(dots$wdata) == FALSE){
        dots$wsigmas <- calcRobustSigmas(dots$wdata,belongmatrix, centers, m = dots$m)
      }
    }


    newbelongmatrix <- update_belongs(data, newcenters, dots)

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
  if(CriterioReached==FALSE){
    warning("The convergence criterion was not reached within the specified number of steps")
  }
  DF <- as.data.frame(newbelongmatrix)
  Groups <- colnames(DF)[max.col(DF, ties.method = "first")]

  results <- c(list(Centers = centers, Belongings = newbelongmatrix,
                  Groups = Groups, Data = data,
                  maxiter = maxiter, tol = tol,
                  isRaster = FALSE),params)

  ## setting the class of the results
  results <- FCMres(results)

  if(dots$noise_cluster){
    results$noise_cluster <- 1 - rowSums(results$Belongings)
    results$delta <- dots$delta
  }

  return(results)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### checking functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Parameter checking function
#'
#' @description Check that the provided parameters are valid
#'
#' @param dots A list of parameters used
#' @param data A numeric and complete dataframe
#' @return A boolean, TRUE if all the tests are passed, FALSE otherwise
#' @importFrom stats complete.cases
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
sanity_check <- function(dots,data){

  ## checking if the parameters have coherent values
  if(is.null(dots$alpha) == FALSE){
    if (dots$alpha < 0){
      stop("alpha parameter must be stricktly superior to 0")
    }
  }

  if(dots$init %in% c("random","kpp") == FALSE){
    stop("the init parameter must be one of 'random' or 'kpp'")
  }

  if(is.null(dots$beta) == FALSE){
    if (dots$beta < 0 | dots$beta > 1){
      stop("beta parameter must be comprised between 0 and 1")
    }
  }

  if(dots$k < 2){
    stop("k must be at least 2")
  }

  if(is.null(dots$seed) == FALSE){
    if (is.integer(dots$seed)){
      stop("the seed parameter must be an integer")
    }
  }

  ## checking if the dataset is complete
  tot <- sum(complete.cases(data)==FALSE)
  if(tot > 0){
    stop("the dataset provided has missing values...")
  }

  ## checking if the dataset is only numeric
  for(col in names(data)){
    if (inherits(data[[col]], c("numeric","integer"))  == FALSE){
      print(paste("the column ",col," is not numeric..."))
      stop("all the columns in the data must be numeric")
    }
  }

  if(!is.null(dots$delta)){
    if(is.na(dots$delta)){
      dots$delta <- NULL
    }
  }

  # checking delta for clustering with noise
  if(dots$noise_cluster & is.null(dots$delta)){
    stop("For clustering with a noise cluster, it is required to set the parameter delta.")
  }
  if(!is.null(dots$delta)){
    if(dots$delta <= 0){
      stop("delta parameter cannot be lower or equal to 0.")
    }
  }

  if(dots$robust & dots$noise_cluster){
    warning("When using the robust version of the algorithm and a noise cluster, consider that the euclidean distance between the observation will be lower because it is normalized. Select the paramter delta accordingly.")
  }

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Adapter functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This set of function is used to provide the parameters in the right order
# to the intermediate functions

#' @title membership matrix calculator for FCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new membership values
#' @keywords internal
belongsFCM <- function(data, centers, dots){
  if(dots$noise_cluster){
    return(calcBelongMatrixNoisy(centers,data,dots$m, dots$delta, dots$sigmas))
  }else{
    return(calcBelongMatrix(centers,data,dots$m, dots$sigmas))
  }

}
#' @title center matrix calculator for FCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param belongmatrix a matrix with the membership values
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new centers
#' @keywords internal
centersFCM <- function(data, centers, belongmatrix, dots){
  return(calcCentroids(data,belongmatrix, dots$m))
}

#' @title membership matrix calculator for GFCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new membership values
#' @keywords internal
belongsGFCM <- function(data,centers,dots){
  if(dots$noise_cluster){
    return(calcFGCMBelongMatrixNoisy(centers,data,dots$m,dots$beta,dots$delta,dots$sigmas))
  }else{
    return(calcFGCMBelongMatrix(centers,data,dots$m,dots$beta,dots$sigmas))
  }

}
#' @title center matrix calculator for GFCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param belongmatrix a matrix with the membership values
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new centers
#' @keywords internal
centersGFCM <- function(data, centers, belongmatrix, dots){
  return(calcCentroids(data,belongmatrix, dots$m))
}

#' @title membership matrix calculator for SFCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new membership values
#' @keywords internal
belongsSFCM <- function(data,centers,dots){
  if(dots$noise_cluster){
    return(calcSFCMBelongMatrixNoisy(centers,data,dots$wdata,dots$m,dots$alpha, dots$delta ,dots$sigmas,dots$wsigmas))
  }else{
    return(calcSFCMBelongMatrix(centers,data,dots$wdata,dots$m,dots$alpha,dots$sigmas,dots$wsigmas))
  }

}
#' @title center matrix calculator for SFCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param belongmatrix a matrix with the membership values
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new centers
#' @keywords internal
centersSFCM <- function(data, centers, belongmatrix, dots){
  return(calcSWFCCentroids(data,dots$wdata,belongmatrix,dots$m,dots$alpha))
}

#' @title membership matrix calculator for SGFCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new membership values
#' @keywords internal
belongsSGFCM <- function(data,centers,dots){
  return(calcSFGCMBelongMatrix(centers,data,dots$wdata,dots$m,dots$alpha,dots$beta,dots$sigmas,dots$wsigmas))
}
#' @title center matrix calculator for SGFCM algorithm
#' @param data a matrix (the dataset used for clustering)
#' @param centers a matrix (the centers of the clusters)
#' @param belongmatrix a matrix with the membership values
#' @param dots a list of other arguments specific to FCM
#' @return a matrix with the new centers
#' @keywords internal
centersSGFCM <- function(data, centers, belongmatrix, dots){
  return(calcSWFCCentroids(data,dots$wdata,belongmatrix,dots$m,dots$alpha))
}
