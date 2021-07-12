#' @title Circular window
#'
#' @description Create a matrix that can be used as a window when working with
#' rasters. It uses a radius to set to 0 the weights of pixels that are farther
#' than this distance. This is helpful to create circular focals.
#'
#' @details The original function come from here: https://scrogster.wordpress.com/2012/10/05/applying-a-circular-moving-window-filter-to-raster-data-in-r/
#' but we reworked it to make it faster and to ensure that the result is matrix with odd dimensions.
#'
#' @param radius The size in meters of the radius of the circular focal
#' @param res The width in meters of a pixel. It is assumed that pixels are squares.
#' @return A binary weight matrix
#' @export
#' @examples
#' # wide of 100 meters for pixels of 2 meters
#' window <- circular_window(100, 2)
#' # setting the center at
#' # row standardisation
#' window_row_std <- window / sum(window)
circular_window<-function(radius, res){
  mat_size <- round(1+(2*radius/res))
  if(mat_size %%2 == 0){
    mat_size <- mat_size + 1
  }
  circ_filter<-matrix(0, nrow=mat_size, ncol=mat_size)
  w <- floor(mat_size / 2)
  dimnames(circ_filter)[[1]]<-seq(-w*res,w*res,by=res)
  dimnames(circ_filter)[[2]]<-seq(-w*res,w*res,by=res)
  m1 <- matrix(rep(as.numeric(dimnames(circ_filter)[[1]]),nrow(circ_filter)),nrow = nrow(circ_filter))
  m2 <- t(m1)
  dists <- sqrt(m1**2 + m2**2)
  window <- ifelse(dists <= radius,1,0)
  return(window)
}
