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
