#' @title Standardizing helper
#'
#' @description Create functions to standardize and unstandardize data
#'
#' @param x a numeric vector or a data.frame with only numeric columns.
#' Non numeric columns are dropped.
#' @return If x was a vector, the function returns a list containing two
#' functions : scale and unscale. The first one is an equivalent of the
#' classical function scale(x, center = TRUE, scale = TRUE). The second
#' can be used to reverse the scaling and get back original units. If x
#' was a data.frame, the same pair of functions is returned inside of
#' a list for each numeric column.
#' @export
#' @examples
#' data(LyonIris)
#' LyonScales <- standardizer(sf::st_drop_geometry(LyonIris))
standardizer <- function(x){

  ## first case, we only have to deal with a numeric vector
  if(inherits(x, "numeric")){
    mu <- mean(x)
    si <- sd(x)
    object <- list(
      "scale" = function(y){return((y-mu)/si)},
      "unscale" = function(y){return(y*si+mu)}
    )
  }else if (inherits(x, "data.frame")){
    test <- sapply(x, is.numeric)
    x <- x[,test]
    object <- lapply(names(x), function(col){
      return(
        list(
          "scale" = function(y){return((y-mean(x[[col]]))/sd(x[[col]]))},
          "unscale" = function(y){return(y*sd(x[[col]])+mean(x[[col]]))}
        )
      )
    })
    names(object) <- names(x)
  }else{
    stop("standardizer can only be used on vector or data.frame")
  }


  return(object)
}
