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
#'         \item Belongings: the final belonging matrix
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
