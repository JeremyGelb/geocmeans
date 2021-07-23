#' @param lag_method A string indicating if a classical lag must be used
#' ("mean") or if a weighted median must be used ("median"). When working with rasters, a function
#' can be given (or a string which will be parsed). It will be applied to all the pixels values
#' in the matrix designated by the parameter window and weighted according to the values of this matrix.
#' Typically, to obtain an average of the pixels in a 3x3 matrix one could use the function sum (or "sum")
#' and set the window as: window <- matrix(1/9,nrow = 3, ncol = 3).
