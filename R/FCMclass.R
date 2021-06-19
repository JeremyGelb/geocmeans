#' @title Instantiate a FCMres object
#'
#' @description Instantiate a FCMres object from a list
#'
#' @param obj A list, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @keywords internal
#' @return A object of class FCMres
#' @examples
#' #This is an internal function, no example provided
FCMres <- function(obj){
  attrs <- names(obj)
  necessary <- c("Centers", "Belongings", "Data")
  if (sum(necessary %in% attrs) < 3){
    stop("The three attributes Centers, Belongings and Data are necessary to create
         an object with class FCMres")
  }

  obj$Centers <- tryCatch(as.matrix(obj$Centers),
           error = function(e)
             print("Obj$Centers must be coercible to a matrix with as.matrix"))

  obj$Belongings <- tryCatch(as.matrix(obj$Belongings),
                      error = function(e)
                        print("Obj$Belongings must be coercible to a matrix with as.matrix"))

  obj$Data <- tryCatch(as.matrix(obj$Data),
               error = function(e)
                 print("Obj$Data must be coercible to a matrix with as.matrix"))


  if("Groups" %in% attrs == FALSE){
    DF <- as.data.frame(obj$newbelongmatrix)
    obj$Groups <- colnames(DF)[max.col(DF, ties.method = "first")]
  }

  class(obj) <- c("FCMres","list")
  return(obj)
}


#' @title is method for FCMres
#'
#' @description Check if an object can me considered as a FCMres object
#'
#' @param x A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @return A boolean, TRUE if x can be considered as a FCMres object, FALSE otherwise
#'   group
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' is.FCMres(result)
is.FCMres <- function(x){
  attrs <- names(x)

  # the object must have this three attributes
  necessary <- c("Centers", "Belongings", "Data")
  if (sum(necessary %in% attrs) < 3){
    return(FALSE)
  }

  # Data, Belongings and Centers must be coercible to matrices
  tryCatch({
    as.matrix(x$Data)
    as.matrix(x$Centers)
    as.matrix(x$Belongings)
    },error = function(e)return(FALSE))

  # the number of rows of Data and Belongings must match
  if(nrow(x$Data) != nrow(x$Belongings)){
    return(FALSE)
  }

  # the number of columns of Data and Centers must match
  if(ncol(x$Data) != ncol(x$Centers)){
    return(FALSE)
  }

  # if we pass all these steps, then the object could be considered as
  # a FCMres
  return(TRUE)

}
