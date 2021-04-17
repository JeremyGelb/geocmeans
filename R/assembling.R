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
      values <- apply(obs,2,function(column){
        return(reldist::wtd.quantile(column,q=0.5,weight=w))
      })
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
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
calcEuclideanDistance <- function(m, v) {
  v <- as.numeric(v)
  alldistances <-colSums((t(m)-v)**2)
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
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
evaluateMatrices <- function(mat1, mat2, tol) {
  mat1 < -as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  differ <- abs(mat1 - mat2)
  diffobs <- apply(differ, 1,max)
  if (length(diffobs[diffobs >= tol]) > 0) {
    return(FALSE)
  } else (return(TRUE))
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

  #checking if the parameters are ok
  sanity_check(dots,data)

  # selection des fonctions de calcul
  if(algo == "FCM"){
    update_belongs <- belongsFCM
    udpdate_centers <- centersFCM

  }else if(algo=="GFCM"){
    update_belongs <- belongsGFCM
    udpdate_centers <- centersGFCM

  }else if (algo=="SFCM"){
    update_belongs <- belongsSFCM
    udpdate_centers <- centersSFCM

  }else if (algo=="SGFCM"){
    update_belongs <- belongsSGFCM
    udpdate_centers <- centersSGFCM
  }

  # selecting the original centers from observations
  if (is.null(dots$seed)==F){
    set.seed(dots$seed)
  }
  centers <- data[sample(nrow(data), k), ]

  # calculating the first belonging matrix
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
  # calculating the most likely group of earch data point
  if(CriterioReached==FALSE){
    warning("The convergence criterion was not reached within the specified number of steps")
  }
  DF <- as.data.frame(newbelongmatrix)
  Groups <- colnames(DF)[max.col(DF, ties.method = "first")]
  return(list(Centers = centers, Belongings = newbelongmatrix,
              Groups = Groups, Data = data))

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

  if(is.null(dots$beta) == FALSE){
    if (dots$beta < 0 | dots$beta > 1){
      stop("beta parameter must be comprised between 0 and 1")
    }
  }

  if(dots$k < 2){
    stop("k must be at leat 2")
  }

  if(is.null(dots$seed) == FALSE){
    if (is.integer(dots$seed)){
      stop("the seed parameter must be an integer")
    }
  }

  ## checking if the dataset is complete
  tot <- sum(complete.cases(data)==F)
  if(tot > 0){
    stop("the dataset provided has missing values...")
  }

  ## checking if the dataset is only numeric
  for(col in names(data)){
    if (class(data[[col]]) %in% c("numeric","integer") == FALSE){
      print(paste("the column ",col," is not numeric..."))
      stop("all the columns in the data must be numeric")
    }
  }

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Adapter functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This set of function is used to provide the parameters in the right order
# to the intermediate functions

belongsFCM <- function(data, centers ,dots){
  return(calcBelongMatrix(centers,data,dots$m))
}
centersFCM <- function(data, centers, belongmatrix, dots){
  return(calcCentroids(data,belongmatrix, dots$m))
}


belongsGFCM <- function(data,centers,dots){
  return(calcFGCMBelongMatrix(centers,data,dots$m,dots$beta))
}
centersGFCM <- function(data, centers, belongmatrix, dots){
  return(calcCentroids(data,belongmatrix, dots$m))
}


belongsSFCM <- function(data,centers,dots){
  return(calcSFCMBelongMatrix(centers,data,dots$wdata,dots$m,dots$alpha))
}
centersSFCM <- function(data, centers, belongmatrix, dots){
  return(calcSWFCCentroids(data,dots$wdata,belongmatrix,dots$m,dots$alpha))
}


belongsSGFCM <- function(data,centers,dots){
  return(calcSFGCMBelongMatrix(centers,data,dots$wdata,dots$m,dots$alpha,dots$beta))
}
centersSGFCM <- function(data, centers, belongmatrix, dots){
  return(calcSWFCCentroids(data,dots$wdata,belongmatrix,dots$m,dots$alpha))
}
