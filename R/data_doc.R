#' @title social and environmental indicators for the Iris of the metropolitan region of Lyon (France)
#'
#' @description A dataset containing social and environmental data for the
#' Iris of Lyon (France)
#' @name LyonIris
#'
#' @format A SpatialPolygonsDataFrame with 506 rows and 32 variables:
#' \describe{
#'   \item{OBJECTID}{a simple OID (integer)}
#'   \item{INSEE_COM}{the code of each commune (factor)}
#'   \item{CODE_IRIS}{the code of each unit area : iris (factor)}
#'   \item{Lden}{the annual daily mean noise exposure values in dB (numeric)}
#'   \item{NO2}{the annual mean of NO2 concentration in ug/m3 (numeric)}
#'   \item{PM25}{the annual mean of PM25 concentration in ug/m3 (numeric)}
#'   \item{PM10}{the annual mean of PM25 concentration in ug/m3 (numeric)}
#'   \item{Pct0_14}{the percentage of people that are 0 to 14 year old (numeric)}
#'   \item{Pct_65}{the percentage of people older than 64 (numeric)}
#'   \item{Pct_Img}{the percentage immigrants (numeric)}
#'   \item{TxChom1564}{the unemployment rate (numeric)}
#'   \item{Pct_brevet}{the percentage of people that obtained the college diploma (numeric)}
#'   \item{NivVieMed}{the median standard of living in euros (numeric)}
#'   \item{VegHautPrt}{the percentage of the iris surface covered by trees (numeric)}
#'   \item{X}{the X coordinate of the center of the Iris (numeric)}
#'   \item{Y}{the Y coordinate of the center of the Iris (numeric)}
#'   ...
#' }
#' @source \url{https://data.grandlyon.com/}
"LyonIris"


#' @title SpatRaster of the bay of Arcachon
#'
#' @description A Landsat 8 image of the bay of Arcachon (France), with a resolution of 30mx30m
#' and 6 bands: blue, green, red, near infrared, shortwave infrared 1 and shortwave infrared 2.
#' The dataset is saved as a Large RasterBrick with the package raster and has the
#' following crs: EPSG:32630. It is provided as a tiff file.
#'
#' @name Arcachon
#'
#' @format A spaRast with 6 bands
#' \describe{
#'   \item{blue}{wavelength: 0.45-0.51}
#'   \item{green}{wavelength: 0.53-0.59}
#'   \item{red}{wavelength: 0.64-0.67}
#'   \item{near infrared}{wavelength: 0.85-0.88}
#'   \item{shortwave infrared}{wavelength: 1.57-1.65}
#'   \item{shortwave infrared}{wavelength: 2.11-2.29}
#' }
#' @source \url{https://earthexplorer.usgs.gov/}
#'
#' @examples
#' # loading directly from file
#' Arcachon <- terra::rast(system.file("extdata/Littoral4_2154.tif", package = "geocmeans"))
#' names(Arcachon) <- c("blue", "green", "red", "infrared", "SWIR1", "SWIR2")
#'
#' # loading with the provided function
#' Arcachon <- load_arcachon()
#' @export
load_arcachon <- function(){
  Arcachon <- terra::rast(system.file("extdata/Littoral4_2154.tif", package = "geocmeans"))
  names(Arcachon) <- c("blue", "green", "red", "infrared", "SWIR1", "SWIR2")
  return(Arcachon)
}


