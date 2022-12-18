context("Testing the behavior of FCMres class")

test_that("Testing the behavior of the FCMres class object : error if not list",{

  expect_error({
    obj <- matrix(0, nrow = , ncol = 3)
    FCMres(obj)
  })

})


test_that("Testing the behavior of the FCMres class object : simple case with vector data",{

  data(LyonIris)
  library(dplyr)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
                     "TxChom1564","Pct_brevet","NivVieMed")
  Data <- sf::st_drop_geometry(LyonIris[AnalysisFields])

  clust <- hclust(dist(Data), method = "ward.D2")
  # getting the groups
  LyonIris$Hclust_groups <- as.character(cutree(clust, k = 4))
  Data$Hclust_groups <- as.character(cutree(clust, k = 4))

  centers <- Data %>%
    group_by(Hclust_groups) %>%
    summarise_all(mean)

  centers <- as.matrix(centers[2:ncol(centers)])
  member_mat <- cat_to_belongings(Data$Hclust_groups)

  Data$Hclust_groups <- NULL

  hclustres <- FCMres(list(
    "Centers" = centers,
    "Belongings" = member_mat,
    "Data" = Data,
    "m" = 1,
    "algo" = "hclust"
  ))

  ## test 1 :
  test1 <- is.FCMres(hclustres)

  ## test 2 :
  Centers <- hclustres$Centers
  hclustres$Centers <- NULL
  test2 <- is.FCMres(hclustres) == FALSE

  ## test 3 :
  err <- expect_error({
    hclustres <- FCMres(list(
    "Centers" = centers,
    "Data" = Data,
    "m" = 1,
    "algo" = "hclust"
  ))
  })

  ## test 4 :
  err <- expect_error({
    hclustres <- FCMres(list(
      "Centers" = centers,
      "Belongings" = "hclust",
      "Data" = Data,
      "m" = 1,
      "algo" = "hclust"
    ))
  })

  ## test 5 :
  err <- expect_error({
    hclustres <- FCMres(list(
      "Centers" = centers,
      "Belongings" = member_mat[1:5,],
      "Data" = Data,
      "m" = 1,
      "algo" = "hclust"
    ))
  })

  ## test 6 :
  member_mat[1,1] <- 5
  err <- expect_warning({
    hclustres <- FCMres(list(
      "Centers" = centers,
      "Belongings" = member_mat,
      "Data" = Data,
      "m" = 1,
      "algo" = "hclust"
    ))
  })

})


test_that("Testing the behavior of the FCMres class object : simple case with raster data",{

  library(terra)
  Arcachon <- terra::rast(system.file("extdata/Littoral4_2154.tif", package = "geocmeans"))
  names(Arcachon) <- c("blue", "green", "red", "infrared", "SWIR1", "SWIR2")

  library(dplyr)
  # loading each raster as a column in a matrix
  # and scale each column
  all_data <- do.call(cbind, lapply(names(Arcachon), function(n){
    rast <- Arcachon[[n]]
    return(terra::values(terra::scale(rast),mat = FALSE))
  }))

  # removing the rows with missing values
  missing <- complete.cases(all_data)
  all_data <- all_data[missing,]

  # applying the kmeans algorithm with 7 groups
  kmean7 <- kmeans(all_data, 7)

  # creating Data (do not forget the standardization)
  Data <- lapply(names(Arcachon), function(n){
    rast <- Arcachon[[n]]
    return(terra::scale(rast))
  })
  names(Data) <- names(Arcachon)

  # creating rasters
  ref_raster <- Arcachon[[1]]

  rasters <- lapply(1:7, function(i){
    # creating a vector with only 0 values
    vals <- rep(0, terra::ncell(ref_raster))
    # filling it with values when the pixels are not NA
    vals[missing] <- ifelse(kmean7$cluster == i,1,0)
    # setting the values in a rasterLayer
    rast <- ref_raster
    terra::values(rast) <- vals
    return(rast)
  })

  # creating centers
  all_data <- as.data.frame(all_data)
  names(all_data) <- names(Arcachon)
  all_data$kmean_groups <- as.character(kmean7$cluster)

  centers <- all_data %>%
    group_by(kmean_groups) %>%
    summarise_all(mean)

  centers <- as.matrix(centers[2:ncol(centers)])

  myFCMres <- FCMres(list(
    "Data" = Data,
    "Centers" = centers,
    "rasters" = rasters,
    "m" = 1,
    "algo" = "kmeans"
  ))

  ## test1 :
  err <- expect_error({
    myFCMres <- FCMres(list(
      "Data" = Data,
      "Centers" = centers,
      "m" = 1,
      "algo" = "kmeans"
    ))
  })

  ## test2 :
  err <- expect_error({
    myFCMres <- FCMres(list(
      "Data" = Data,
      "Centers" = centers,
      "rasters" = centers,
      "m" = 1,
      "algo" = "kmeans"
    ))
  })

})

