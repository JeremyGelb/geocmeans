#' @param window If data is a list of rasters, then a window must be specified instead of
#' a list.w object. It will be used to calculate a focal function on each raster. The
#' window must be a square numeric matrix with odd dimensions (such 3x3). The values in
#' the matrix indicate the weight to give to each pixel and the centre of the matrix is
#' the centre of the focal function.
