
#' @title Global Moran I for raster
#'
#' @description Calculate the global Moran I for a numeric raster
#'
#' @param rast A RasterLayer or a matrix
#' @param window The window defining the neighbour weights
#' @return A float: the global Moran I
#' @examples
#' # this is an internal function, no example provided
calc_moran_raster <- function(rast, window){
  if(class(rast)[[1]] == "RasterLayer"){
    mat <- raster::as.matrix(rast)
  }else if (class(rast)[[1]] == "matrix"){
    mat <- rast
  }else{
    stop("rast parameter must be on of matrix or RasterLayer")
  }
  if(class(window)[[1]] == "matrix"){
    fun <- moranI_matrix_window
  }else if (class(window)[[1]] == "numeric"){
    fun <- moranI_matrix
  }else{
    stop("window parameter must be an integer or a matrix")
  }
  val <- fun(mat, window)
  return(val)
}


#' @title Local Moran I for raster
#'
#' @description Calculate the Local Moran I for a numeric raster
#'
#' @param rast A RasterLayer or a matrix
#' @param window The window defining the neighbour weights
#' @return A vector: the local Moran I
#' @examples
#' # this is an internal function, no example provided
calc_local_moran_raster <- function(rast, window){
  if(class(rast)[[1]] == "RasterLayer"){
    mat <- raster::as.matrix(rast)
  }else if (class(rast)[[1]] == "matrix"){
    mat <- rast
  }else{
    stop("rast parameter must be on of matrix or RasterLayer")
  }
  if(class(window)[[1]] == "matrix"){
    window <- window
  }else if (class(window)[[1]] == "numeric"){
    w <- 1+2*window
    window <- matrix(1, nrow = w, ncol = w)
    window[ceiling(w/2),ceiling(w/2)] <- 0
  }else{
    stop("window parameter must be an integer or a matrix")
  }
  vals <- local_moranI_matrix_window(mat, window)

  if(class(rast)[[1]] == "RasterLayer"){
    raster::values(rast) <- vals
    return(rast)
  }else{
    return(matrix(vals, ncol = ncol(mat), nrow = nrow(mat)))
  }

}
