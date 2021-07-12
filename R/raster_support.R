# # https://rdrr.io/bioc/OLIN/man/ma.matrix.html
# # https://www.rdocumentation.org/packages/raster/versions/3.4-10/topics/focal
# # example de jeu de donnees : https://www.r-exercises.com/2018/02/28/advanced-techniques-with-raster-data-part-1-unsupervised-classification/?__cf_chl_jschl_tk__=703ae8c7214fe36addfe12006afebc7442de2077-1624418523-0-AQOrJJeKBwWikYZcyrzGFdsNhq2ERFu4EIyPoPuvPuVd77RyG-ppxAcA-KJizHAhNWTVXc7crst9Hs-IBtbOl1yfjSzVGcRFzlsof2sCbYRPk4F5GrKrhn64zm5fOGhiMi25RcoqrnuYbjLpu5QFzyVtBlWiSvQ_IF1SO27b388v_IXyW5G-BHA1htk5vJeHBHnUZDDmDQNUeVhRDzgVrZVjCSP1akW5P854dj0ILOKNcYoqrAQVxvpUaonarM38I6QiQaTcgTpnXxSPXkNYaOAmLkReEjcpYXTTT4gETRH5Ymnd8fDElBBq9pkmvKxjGc-Sx5-Y15o-iAGKenGYWQSdXrypzaML-mMAJSSM7GEgL2a6dJ1HRtOHcxwmVAW6MuUPR2cabYkEzyTlWgfZsWk_l92jaiaqo_7KCo92WZRxWk3s3f3CfKhLenCglEjXJ-F9Rkn-LDPD-Iuh4NDABe8YHHs_F7Gbrb9TL__ii50bNlx9rN03RBV4H4ZZcRaSiw

# # # # # #TODO : add waiter package to have a nice waiting pannel !
# library(raster)
# raster1 <- raster(".Rproj.user/raster_data_test/estuary_LC.tif")
#
# data <- lapply(1:nbands(raster1), function(i){
#   myband <- raster(".Rproj.user/raster_data_test/estuary_LC.tif", band = i)
#   return(myband)
# })
#
# #names(data) <- paste("band",1:nbands(raster1), sep="")
#
# window <- matrix(1, nrow = 3, ncol = 3)
# window <- window / sum(window)
#
# lag_method  <- "sum"
# standardize <- TRUE
# na.rm <- TRUE
#
# k = 5
# m = 1.5
# alpha = 1
# lag_method=sum
# maxiter = 1000
# tol = 0.001
# standardize = TRUE
# verbose = TRUE
# seed = 123
# init = "random"
#
# object <- SFCMeans(data = data, k = k, m = m,
#                    window = window,
#                    alpha = alpha, lag_method = lag_method, maxiter = maxiter,
#                    tol = tol, standardize = standardize)

# titi <- select_parameters("SFCM", data = data, k = k, m = m, alpha = c(1,2),
#                           lag_method = "sum", maxiter = maxiter, tol = tol,
#                           standardize = standardize, window  = window
#                           )


#' @title Check the shape of a window
#'
#' @description Check is a window is squarred and have odd dimensions
#'
#' @param w A matrix
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
check_window <- function(w){
  dims <- dim(w)
  if(dims[[1]] != dims[[2]]){
    stop("The window provided must be a square matrix")
  }
  if(dims[[1]] %%2 == 0){
    stop("The width and height of the window must be an odd number")
  }
}

#' @title Raster data preparation
#'
#' @description Prepare a raster dataset
#'
#' @param dataset A list of rasters
#' @param w The window to use in the focal function
#' @param fun the function to use as the focal function
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduced (default = True)
#' @return A list with the required elements to perform clustering
#' @importFrom raster focal
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
input_raster_data <- function(dataset, w = NULL, fun = sum, standardize = TRUE){

  if(is.null(w) == FALSE){

    check_window(w)

    if(class(fun) != "function"){

      tryCatch(fun <- eval(parse(text=fun)),
               error = function(e)
                 print("When using rasters, the parameter lag_method must be a function or a string
                 that can be parsed to a function like sum, mean, min, max ...,
                 Note that this function is applied to the pixels values multiplied by the weights in the window.")
               )
    }
  }

  refdim <- dim(dataset[[1]])
  for(band in dataset){
    dim2 <- dim(band)
    if(sum(abs(dim2 - refdim)) != 0){
      stop("The rasters provided must have EXACTLY the same dimensions.")
    }
  }

  if(is.null(names(dataset))){
    okname <- paste("band",1:length(dataset), sep="")
  }else{
    okname <- names(dataset)
  }

  if(standardize){
    for(i in 1:length(dataset)){
      dataset[[i]] <- raster::scale(dataset[[i]])
    }
  }

  #step1 : prepare the regular dataset
  datamatrix <- do.call(cbind,lapply(dataset, function(band){
    raster::values(band)
  }))

  #step2 : only keep the non-missing valuez
  missing_pxl <- complete.cases(datamatrix)
  data_class <- datamatrix[missing_pxl,]
  colnames(data_class) <- okname

  #step3 : calculating Wdata if necessary
  if(is.null(w) == FALSE){
    Wdatamatrix <- do.call(cbind,lapply(dataset, function(band){
      wraster <- focal(band, w, sum, na.rm = TRUE, pad = TRUE)
      return(raster::values(wraster))
    }))

    Wdata_class <- Wdatamatrix[missing_pxl,]
    colnames(Wdata_class) <- okname

  }else{
    Wdata_class <- NULL
  }

  return(list("data" = data_class,
              "wdata" = Wdata_class,
              "missing" = missing_pxl,
              "rst" = dataset[[1]],
              "window" = w
              ))
}



#' @title Raster result transformation
#'
#' @description Adapt the results if a raster is used
#'
#' @param object A FCMres object
#' @param missing A boolean indicating which pixels have no missing values
#' @param rst A raster object used as template to structure the results
#' @return A FCMres object with isRaster = TRUE
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
output_raster_data <- function(object, missing, rst){
  # object is created by a FCM like function
  # missing is indicating which pixels are complete
  # rst is a basic raster to duplicate it

  #step1 : creating a raster for each cluster
  rasters <- lapply(1:ncol(object$Belongings), function(i){
    rst2 <- rst
    vec <- rep(NA,times = ncell(rst))
    vec[missing] <- object$Belongings[,i]
    raster::values(rst2) <- vec
    return(rst2)
  })

  names(rasters) <- paste("group", 1:ncol(object$Belongings), sep = "")

  #step2 : adding the most likely group
  rst2 <- rst
  vec <- rep(NA,times = ncell(rst))
  DF <- as.data.frame(object$Belongings)
  vec[missing] <- max.col(DF, ties.method = "first")
  raster::values(rst2) <- vec
  rasters[["Groups"]] <- rst2

  object$isRaster <- TRUE
  object$rasters <- rasters
  object$missing <- missing
  return(object)
}

