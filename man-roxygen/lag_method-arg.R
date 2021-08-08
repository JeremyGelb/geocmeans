#' @param lag_method A string indicating if a classical lag must be used
#' ("mean") or if a weighted median must be used ("median"). When working with rasters, a function
#' can be given (or a string which will be parsed). It will be applied to all the pixels values
#' in the matrix designated by the parameter window and weighted according to the values of this matrix.
#' Typically, to obtain an average of the pixels in a 3x3 matrix one could use the function sum (or "sum")
#' and set the window as: window <- matrix(1/9,nrow = 3, ncol = 3). There is one special case when working 
#' with rasters: one can specify "nl" (standing for non-local) which calculated a lagged version of the 
#' input rasters, using the inverse of the euclidean distance as spatial weights (see the section Advanced 
#' examples in the vignette introduction for more details).
