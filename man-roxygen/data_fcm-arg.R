#' @param data A dataframe with only numerical variables. Can also be a list of
#' rasters (produced by the package raster). In that case, each raster is
#' considered as a variable and each pixel is an observation. Pixels with NA
#' values are not used during the classification.