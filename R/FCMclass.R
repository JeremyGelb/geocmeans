#' @title Instantiate a FCMres object
#'
#' @description Instantiate a FCMres object from a list
#'
#' @details Creating manually a FCMres object can be handy to use geocmeans functions
#' on results from external algorithms. The list given to function FCMres must
#' contain 5 necessary parameters:
#' \itemize{
#'         \item Centers: a dataframe or matrix describing the final centers of the groups
#'         \item Belongings: a membership matrix
#'         \item Data: the dataset used to perform the clustering.
#' It must be a dataframe or a matrix. If a list is given, then the function
#' assumes that the classification occured on rasters (see information below)
#'         \item m: the fuzyness degree (1 if hard clustering is used)
#'         \item algo: the name of the algorithm used
#' }
#'
#' Note that the S3 method predict is available only for object created with the functions
#' CMeans, GCMeans, SFCMeans, SGFCMeans.
#'
#' When working with rasters, Data must be a list of rasters, and a second list of rasters
#' with the membership values must be provided is en extra slot named "rasters". In that
#' case, Belongings has not to be defined and will be created automatically.
#'
#' Warning: the order of the elements is very important. The first row in the
#' matrix "Centers", and the first column in the matrix "Belongings" must both
#' be related to the same group and so on. When working with raster data, the
#' first row in the matrix "Centers" must also match with the first rasterLayer
#' in the list "rasters".
#'
#' @param obj A list, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @keywords internal
#' @return An object of class FCMres
#' @export
#' @examples
#' #This is an internal function, no example provided
FCMres <- function(obj){

  #if(class(obj)!="list"){
  if (inherits(obj,"list") == FALSE){
    stop("A list is required to create a FCMres object.")
  }

  #rasterMode <- class(obj$Data)[[1]] == "list"
  rasterMode <- inherits(obj$Data, c("list"))
  if(rasterMode){
    necessary <- c("Centers", "Data", "m", "algo", "rasters")
    attrs <- names(obj)
    if (sum(necessary %in% attrs) < 5){
      stop("The attributes Centers, Data, m, algo and rasters are necessary to create
         an object with class FCMres from raster data (obj$Data is a list, see
         details in help(FCMres))")
    }
  }else{
    necessary <- c("Centers", "Belongings", "Data", "m", "algo")
    attrs <- names(obj)
    if (sum(necessary %in% attrs) < 5){
      stop("The attributes Centers, Data, m, algo and Belongings are necessary to create
         an object with class FCMres")
    }
    obj$Belongings <- tryCatch(as.matrix(obj$Belongings),
                               error = function(e)
                                 print("Obj$Belongings must be coercible to a matrix with as.matrix"))
    if(nrow(obj$Belongings) != nrow(obj$Data)){
      stop("The number of rows in obj$Data and obj$Belongings must be the same.")
    }
  }

  obj$Centers <- tryCatch(as.matrix(obj$Centers),
           error = function(e)
             print("Obj$Centers must be coercible to a matrix with as.matrix"))



  if(rasterMode){
    # if we are in raster mode, the data are provided as a list of rasterLayers
    check_raters_dims(obj$Data)

    #step1 : prepare the regular dataset
    datamatrix <- do.call(cbind,lapply(obj$Data, function(band){
      terra::values(band, mat = FALSE)
    }))

    #step2 : only keep the non-missing values
    missing_pxl <- complete.cases(datamatrix)
    data_class <- datamatrix[missing_pxl,]
    obj$Data <- data_class
    obj$missing <- missing_pxl
    obj$isRaster <- TRUE
    if(inherits(obj$rasters, "list") == FALSE){
    #if(is.null(obj$rasters) | class(obj$rasters)!="list"){
      stop('When using raster data, the slot "raster" in obj must contains a list of the
           rasterLayers representing the membership values')
    }
    #step3 : set up the membership matrix
    belongmat <- do.call(cbind,lapply(obj$rasters, function(band){
      terra::values(band, mat = FALSE)
    }))
    belongmat <- belongmat[missing_pxl,]
    obj$Belongings <- belongmat

    #step4 pimp a bit the raster object
    old_raster <- obj$rasters
    names(old_raster) <- paste("group",1:length(old_raster))

    rst2 <- old_raster[[1]]
    vec <- rep(NA,times = terra::ncell(rst2))
    DF <- as.data.frame(obj$Belongings)
    vec[obj$missing] <- max.col(DF, ties.method = "first")
    terra::values(rst2) <- vec
    old_raster[["Groups"]] <- rst2
    obj$rasters <- old_raster

  }else{
    obj$Data <- tryCatch(as.matrix(obj$Data),
                         error = function(e)
                           print("Obj$Data must be coercible to a matrix with as.matrix"))
    obj$isRaster <- FALSE
  }

  obj$k <- ncol(obj$Belongings)


  if("Groups" %in% attrs == FALSE){
    DF <- as.data.frame(obj$Belongings)
    obj$Groups <- colnames(DF)[max.col(DF, ties.method = "first")]
  }

  # A quick check of the membership matrix
  test1 <- any(round(rowSums(obj$Belongings),8) != 1)
  test2 <- is.null(obj$noise_cluster)
  if(test1 & test2){
    warning("some rows in the membership matrix does not sum up to 1... This should be checked")
  }

  class(obj) <- c("FCMres","list")
  return(obj)
}


#' @title is method for FCMres
#'
#' @description Check if an object can be considered as a FCMres object
#'
#' @param object A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @param class2 Character string giving the names of the classe to test (usually "FCMres")
#' @return A boolean, TRUE if x can be considered as a FCMres object, FALSE otherwise
#'   group
#' @importFrom methods is
#' @method is FCMres
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' is(result, "FCMres")
is.FCMres <- function(object, class2 = 'FCMres'){
  attrs <- names(object)

  # the object must have this three attributes
  necessary <- c("Centers", "Belongings", "Data", "m", "algo")
  if (sum(necessary %in% attrs) < 5){
    return(FALSE)
  }

  # Data, Belongings and Centers must be coercible to matrices
  tryCatch({
    as.matrix(object$Data)
    as.matrix(object$Centers)
    as.matrix(object$Belongings)
    },error = function(e)return(FALSE))

  # the number of rows of Data and Belongings must match
  if(nrow(object$Data) != nrow(object$Belongings)){
    return(FALSE)
  }

  # the number of columns of Data and Centers must match
  if(ncol(object$Data) != ncol(object$Centers)){
    return(FALSE)
  }

  # if we pass all these steps, then the object could be considered as
  # a FCMres
  return(TRUE)
}


#' @title print method for FCMres
#'
#' @description print a FCMres object
#'
#' @param x A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @param ... not used
#' @return A boolean, TRUE if x can be considered as a FCMres object, FALSE otherwise
#'   group
#' @method print FCMres
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' result <- CMeans(dataset, k = 5, m = 1.5, standardize = TRUE)
#' print(result, "FCMres")
print.FCMres <- function(x, ...){ #nocov start
  if(x$isRaster){
    part1 <- "FCMres object constructed from a raster dataset"
    part2 <- paste("dimension of the raster : ",x$rasters[[1]]@nrows,"x",x$rasters[[1]]@ncols,
                   " and ", nrow(x$Data), " variables",sep="")
  }else{
    part1 <- "FCMres object constructed from a dataframe"
    part2 <- paste(nrow(x$Data), " observations and ", ncol(x$Data), " variables", sep="")
  }
  elements <- c("algo","k", "m", "alpha", "beta")
  params <- paste(paste(elements, x[elements], sep = "="), collapse = " ; ")
  part3 <- paste("parameters of the classification: ",params, sep = "")
  cat(part1,part2,part3,sep="\n")
  invisible(x)
} #nocov end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####                      SUMMARY S3 METHOD                          ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Descriptive statistics by group
#'
#' @description Calculate some descriptive statistics of each group
#'
#' @param data The original dataframe used for the classification
#' @param belongmatrix A membership matrix
#' @param weighted A boolean indicating if the summary statistics must use the
#'   membership matrix columns as weights (TRUE) or simply assign each
#'   observation to its most likely cluster and compute the statistics on each
#'   subset (FALSE)
#' @param dec An integer indicating the number of digits to keep when rounding
#' (default is 3)
#' @param silent A boolean indicating if the results must be printed or silently returned
#' @return A list of length k (the number of group). Each element of the list is
#'   a dataframe with summary statistics for the variables of data for each
#'   group
#' @export
#' @importFrom dplyr %>%
#' @importFrom grDevices rgb
#' @importFrom stats quantile sd weighted.mean
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' summarizeClusters(dataset, result$Belongings)
summarizeClusters <- function(data, belongmatrix, weighted = TRUE, dec = 3, silent=TRUE) { #nocov start
  belongmatrix <- as.data.frame(belongmatrix)
  if (weighted) {
    Summaries <- lapply(1:ncol(belongmatrix), function(c) {
      W <- belongmatrix[, c]
      Values <- apply(data, 2, function(x) {
        Q5 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.05, na.rm = TRUE, weight = W),dec))
        Q10 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.1, na.rm = TRUE, weight = W),dec))
        Q25 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.25, na.rm = TRUE, weight = W),dec))
        Q50 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.5, na.rm = TRUE, weight = W),dec))
        Q75 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.75, na.rm = TRUE, weight = W),dec))
        Q90 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.9, na.rm = TRUE, weight = W),dec))
        Q95 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.95, na.rm = TRUE, weight = W),dec))
        Mean <- as.numeric(round(weighted.mean(x, W),dec))
        Std <- round(sqrt(reldist::wtd.var(x, weight=W)),dec)
        N <-
          return(list(Q5 = Q5, Q10 = Q10, Q25 = Q25, Q50 = Q50,
                      Q75 = Q75, Q90 = Q90, Q95 = Q95,
                      Mean = Mean, Std = Std))
      })
      DF <- do.call(cbind, Values)
      if(silent==FALSE){
        print(paste("Statistic summary for cluster ", c, sep = ""))
        print(DF)
      }

      return(DF)
    })
    names(Summaries) <- paste("Cluster_", c(1:ncol(belongmatrix)), sep = "")
    return(Summaries)

  } else {
    Groups <- colnames(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    data$Groups <- Groups
    Summaries <- lapply(unique(data$Groups), function(c) {
      DF <- subset(data, data$Groups == c)
      DF$Groups <- NULL
      Values <- apply(DF, 2, function(x) {
        Q5 <- as.numeric(round(quantile(x, probs = 0.05, na.rm = TRUE),dec))
        Q10 <- as.numeric(round(quantile(x, probs = 0.1, na.rm = TRUE),dec))
        Q25 <- as.numeric(round(quantile(x, probs = 0.25, na.rm = TRUE),dec))
        Q50 <- as.numeric(round(quantile(x, probs = 0.5, na.rm = TRUE),dec))
        Q75 <- as.numeric(round(quantile(x, probs = 0.75, na.rm = TRUE),dec))
        Q90 <- as.numeric(round(quantile(x, probs = 0.9, na.rm = TRUE),dec))
        Q95 <- as.numeric(round(quantile(x, probs = 0.95, na.rm = TRUE),dec))
        Mean <- round(mean(x),dec)
        Std <- round(sd(x),dec)
        return(list(Q5 = Q5, Q10 = Q10, Q25 = Q25, Q50 = Q50,
                    Q75 = Q75, Q90 = Q90, Q95 = Q95,
                    Mean = Mean, Std = Std))
      })
      DF <- do.call(cbind, Values)
      if(silent==FALSE){
        print(paste("Statistic summary for cluster ", c, sep = ""))
        print(DF)
      }

      return(DF)
    })
    names(Summaries) <- paste("Cluster_", c(1:ncol(belongmatrix)), sep = "")
    return(Summaries)
  }
}


#' @title Summary method for FCMres
#'
#' @description Calculate some descriptive statistics of each group of a
#' FCMres object
#'
#' @param object A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @param data A dataframe to use for the summary statistics instead of obj$data
#' @param weighted A boolean indicating if the summary statistics must use the
#'   membership matrix columns as weights (TRUE) or simply assign each
#'   observation to its most likely cluster and compute the statistics on each
#'   subset (FALSE)
#' @param dec An integer indicating the number of digits to keep when rounding
#' (default is 3)
#' @param silent A boolean indicating if the results must be printed or silently returned
#' @param ... Not used
#' @return A list of length k (the number of group). Each element of the list is
#'   a dataframe with summary statistics for the variables of data for each
#'   group
#' @method summary FCMres
#' @export
#' @importFrom dplyr %>%
#' @importFrom grDevices rgb
#' @importFrom stats quantile sd weighted.mean
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' summary(result)
summary.FCMres <- function(object, data = NULL, weighted = TRUE, dec = 3, silent=TRUE, ...) {
  if(is.null(data)){
    df <- object$Data
  }else{
    df <- as.matrix(data)
  }
  return(summarizeClusters(df, object$Belongings, weighted, dec, silent))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####                      PREDICT S3 METHOD                          ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Predict matrix membership for new observations
#'
#' @description Function to predict the membership matrix of a new set of observations
#'
#' @template FCMresobj-arg
#' @param new_data A DataFrame with the new observations or a list of rasters if
#' object$isRaster is TRUE
#' @template nblistw-arg
#' @template window-arg
#' @param standardize  A boolean to specify if the variable must be centered and
#'   reduced (default = True)
#' @param ... not used
#'
#' @return A numeric matrix with the membership values for each new observation. If
#' rasters were used, return a list of rasters with the membership values.
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#'
#' # rescaling all the variables used in the analysis
#' for (field in AnalysisFields) {
#'     LyonIris[[field]] <- scale(LyonIris[[field]])
#' }
#'
#' # doing the initial clustering
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = FALSE)
#'
#' # using a subset of the original dataframe as "new data"
#' new_data <- LyonIris[c(1, 27, 36, 44, 73),]
#' new_dataset <- sf::st_drop_geometry(new_data[AnalysisFields])
#' new_nb <- spdep::poly2nb(new_data,queen=TRUE)
#' new_Wqueen <- spdep::nb2listw(new_nb,style="W")
#'
#' # doing the prediction
#' predictions <- predict_membership(result, new_dataset, new_Wqueen, standardize = FALSE)
predict_membership <- function(object, new_data, nblistw = NULL, window = NULL, standardize = TRUE, ...){

  if(object$algo %in% c("FCM", "GFCM", "SFCM", "SGFCM") == FALSE){
    stop('prediction can only be performed for FCMres object
         if algo is one of "FCM", "GFCM", "SFCM", "SGFCM"')
  }


  results <- object
  # if we are in raster mode, we need to do some conversions.
  if(results$isRaster == TRUE){

    # let us check if all the rasters have the same dims
    check_raters_dims(new_data)

    old_data <- new_data
    new_data_tot <- do.call(cbind,lapply(new_data, function(rast){
      return(terra::values(rast, mat = FALSE))
    }))
    missing <- complete.cases(new_data_tot)
    new_data <- new_data_tot[missing,]
  }
  ## standardize the data if required
  if (standardize){
    for (i in 1:ncol(new_data)) {
      new_data[, i] <- scale(new_data[, i])
    }
  }

  ## calculating the lagged dataset if required
  if(results$algo %in% c("SFCM", "SGFCM")){
    if(results$isRaster == FALSE){
      if(is.null(nblistw)){
        stop("With a spatial clustering method like SFCM or SGFCM, the spatial matrix listw
                 associated with the new dataset must be provided")
      }
      wdata <- calcLaggedData(new_data, nblistw, results$lag_method)
    }else{
      if(is.null(window)){
        stop("With a spatial clustering method like SFCM or SGFCM, the spatial matrix window
                 to use on the new dataset must be provided")
      }
      fun <- results$lag_method
      #if(class(fun) != "function"){
      if(inherits(fun, "function")==FALSE){

        tryCatch(fun <- eval(parse(text=fun)),
                 error = function(e)
                   print("When using rasters, the parameter lag_method must be a function or a string
                 that can be parsed to a function like sum, mean, min, max ...,
                 Note that this function is applied to the pixels values multiplied by the weights in the window.")
        )
      }

      wdata_total <- do.call(cbind,lapply(dataset, function(band){
        wraster <- terra::focal(band, window, fun, na.rm = TRUE, pad = TRUE)
        return(terra::values(wraster, mat = FALSE))
      }))
      wdata <- wdata_total[missing,]
    }

  }
  noise_cluster <- is.numeric(results$noise_cluster)
  if(results$robust == FALSE){
    sigmas <- rep(1, nrow(results$Centers))
    wsigmas <- rep(1, nrow(results$Centers))
  }else{
    sigmas <- calcRobustSigmas(data = as.matrix(new_data),
                               belongmatrix = results$Belongings,
                               centers = as.matrix(results$Centers),
                               m = results$m
                               )
    if(results$algo %in% c("SFCM", "SGFCM")){
      print("hello there !")
      wsigmas <- calcRobustSigmas(data = as.matrix(wdata),
                                  belongmatrix = results$Belongings,
                                  centers = as.matrix(results$Centers),
                                  m = results$m
      )
    }
  }

  ## selecting the appropriate function for prediction
  if(results$algo == "FCM"){
    if(!noise_cluster){
      pred_values <- calcBelongMatrix(as.matrix(results$Centers), as.matrix(new_data),
                                      m = results$m, sigmas = sigmas)
    }else{
      pred_values <- calcBelongMatrixNoisy(as.matrix(results$Centers), as.matrix(new_data),
                                      m = results$m, sigmas = sigmas,
                                      delta = results$delta)
    }

  }else if(results$algo == "GFCM"){
    if(!noise_cluster){
      pred_values <- calcFGCMBelongMatrix(as.matrix(results$Centers),
                                          as.matrix(new_data),
                                          m = results$m, beta = results$beta,
                                          sigmas = sigmas)
    }else{
      pred_values <- calcFGCMBelongMatrixNoisy(as.matrix(results$Centers),
                                          as.matrix(new_data),
                                          m = results$m, beta = results$beta,
                                          delta = results$delta,
                                          sigmas = sigmas)
    }

  }else if(results$algo == "SFCM"){
    if(!noise_cluster){
      pred_values <- calcSFCMBelongMatrix(as.matrix(results$Centers),
                                          as.matrix(new_data),
                                          wdata = as.matrix(wdata),
                                          m = results$m,
                                          alpha = results$alpha,
                                          wsigmas = wsigmas,
                                          sigmas = sigmas)
    }else{
      pred_values <- calcSFCMBelongMatrixNoisy(as.matrix(results$Centers),
                                          as.matrix(new_data),
                                          wdata = as.matrix(wdata),
                                          m = results$m, alpha = results$alpha,
                                          delta = results$delta,
                                          wsigmas = wsigmas,
                                          sigmas = sigmas)
    }

  }else if(results$algo == "SGFCM"){
    if(!noise_cluster){
      pred_values <- calcSFGCMBelongMatrix(as.matrix(results$Centers),
                                           as.matrix(new_data),
                                           wdata = as.matrix(wdata),
                                           m = results$m, alpha = results$alpha,
                                           beta = results$beta,
                                           wsigmas = wsigmas,
                                           sigmas = sigmas)
    }else{
      pred_values <- calcSFGCMBelongMatrixNoisy(as.matrix(results$Centers),
                                           as.matrix(new_data),
                                           wdata = as.matrix(wdata),
                                           m = results$m, alpha = results$alpha,
                                           beta = results$beta,
                                           sigmas = sigmas,
                                           wsigmas = wsigmas,
                                           delta = results$delta)
    }

  }

  if(results$isRaster == FALSE){
    return(pred_values)
  }else{
    nc <- terra::ncell(new_data[[1]])
    rasters_membership <- lapply(1:ncol(pred_values), function(i){
      rast <- old_data[[1]]
      vec <- rep(NA,times = nc)
      vec[missing] <- pred_values[,i]
      terra::values(rast, mat = FALSE) <- vec
      return(rast)
    })
    names(rasters_membership) <- paste("group",1:ncol(pred_values), sep = "")
    return(rasters_membership)
  }

}


#' @title Predict method for FCMres object
#' @description Function to predict the membership matrix of a new set of observations
#' @template FCMresobj-arg
#' @param new_data A DataFrame with the new observations
#' @template nblistw-arg
#' @template window-arg
#' @param standardize  A boolean to specify if the variable must be centred and
#'   reduced (default = True)
#' @param ... not used
#' @return A numeric matrix with the membership values for each new observation
#' @method predict FCMres
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#'
#' # rescaling all the variables used in the analysis
#' for (field in AnalysisFields) {
#'     LyonIris[[field]] <- scale(LyonIris[[field]])
#' }
#'
#' # doing the initial clustering
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = FALSE)
#'
#' # using a subset of the original dataframe as "new data"
#' new_data <- LyonIris[c(1, 27, 36, 44, 73),]
#' new_dataset <- sf::st_drop_geometry(new_data[AnalysisFields])
#' new_nb <- spdep::poly2nb(new_data,queen=TRUE)
#' new_Wqueen <- spdep::nb2listw(new_nb,style="W")
#'
#' # doing the prediction
#' predictions <- predict(result, new_dataset, new_Wqueen, standardize = FALSE)
predict.FCMres <- function(object, new_data, nblistw = NULL, window = NULL, standardize = TRUE, ...){
  return(predict_membership(object, new_data, nblistw, window, standardize, ...))
}



#' @title Plot method for FCMres object
#' @description Method to plot the results of a FCM.res object
#'
#' @details This S3 method is a simple dispatcher for the functions barPlots,
#' violinPlots and spiderPlots. To be able to use all their specific
#' parameters, one can use them directly.
#'
#' @param x A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans. Can also be a simple membership matrix.
#' @param type A string indicating the type of plot to show. Can be one of
#' "bar", "violin", or "spider". Default is spider.
#' @param ... not used
#' @return a ggplot2 object, a list, or NULL, depending on the type of plot requested
#' @method plot FCMres
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#'
#' # rescaling all the variables used in the analysis
#' for (field in AnalysisFields) {
#'     LyonIris[[field]] <- scale(LyonIris[[field]])
#' }
#'
#' # doing the initial clustering
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = FALSE)
#'
#' plot(result, type = "spider")
plot.FCMres <- function(x, type = "spider", ...){
  if(type %in% c("bar", "violin","spider") == FALSE){
    stop('the parameter type must be one of "bar", "violin","spider" to plot a FCM.res object')
  }
  if(type == "bar"){
    return(barPlots(x$Data, x$Belongings))
  }else if(type == "violin"){
    return(violinPlots(x$Data,x$Groups))
  }else if(type == "spider"){
    return(spiderPlots(x$Data,x$Belongings))
  }
}#nocov end
