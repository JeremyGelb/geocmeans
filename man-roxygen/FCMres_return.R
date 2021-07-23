#' @return An S3 object of class FCMres with the following slots
#' \itemize{
#'         \item Centers: a dataframe describing the final centers of the groups
#'         \item Belongings: the final membership matrix
#'         \item Groups: a vector with the names of the most likely group for each observation
#'         \item Data: the dataset used to perform the clustering (might be standardized)
#'         \item isRaster: TRUE if rasters were used as input data, FALSE otherwise
#'         \item k: the number of groups
#'         \item m: the fuzyness degree
#'         \item alpha: the spatial weighting parameter (if SFCM or SGFCM)
#'         \item beta: beta parameter for generalized version of FCM (GFCM or SGFCM)
#'         \item algo: the name of the algorithm used
#'         \item rasters: a list of rasters with membership values and the most likely group (if rasters were used)
#'         \item missing: a boolean vector indicating raster cell with data (TRUE) and with NA (FALSE) (if rasters were used)
#'         \item maxiter: the maximum number of iterations used
#'         \item tol: the convergence criterio
#'         \item lag_method: the lag function used (if SFCM or SGFCM)
#'         \item nblistw: the neighbours list used (if vector data were used for SFCM or SGFCM)
#'         \item window: the window used (if raster data were used for SFCM or SGFCM)
#' }
#' @export
