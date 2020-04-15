context("belonging matrices")

test_that("SFCM belonging matrices must have a rowSums equal to 1",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")
  result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
  x <- round(1-rowSums(result$Belongings),6)
  expect_equal(sum(x),0)
})

test_that("SFCM belonging matrices must have a rowSums equal to 1, with high alpha",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")
  result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 8, standardize = TRUE)
  x <- round(1-rowSums(result$Belongings),6)
  expect_equal(sum(x),0)
})

test_that("SFCM belonging matrices must have a rowSums equal to 1, with 0 alpha",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")
  result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 0, standardize = TRUE)
  x <- round(1-rowSums(result$Belongings),6)
  expect_equal(sum(x),0)
})


test_that("When alpha is 0, the belonging matrix should be the same for the two methods",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  data <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")

  for (i in 1:ncol(data)) {
    data[, i] <- scale(data[, i])
  }

  wdata <- data
  for (Name in names(data)) {
    wdata[[Name]] <- spdep::lag.listw(Wqueen, data[[Name]])
  }

  centers <- data[sample(nrow(data), 4), ]

  belongMat1 <- calcBelongMatrix(centers, data, 1.5)
  belongMat2 <- calcSFCMBelongMatrix(centers, data, wdata, 1.5,0)

  diff1 <- belongMat1 - belongMat2
  tot1 <- round(sum(diff1),8)
  expect_equal(tot1,0)
})


test_that("When alpha is 0, the centers should be the same for the two methods",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  data <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")

  for (i in 1:ncol(data)) {
    data[, i] <- scale(data[, i])
  }

  wdata <- data
  for (Name in names(data)) {
    wdata[[Name]] <- spdep::lag.listw(Wqueen, data[[Name]])
  }

  centers <- data[sample(nrow(data), 4), ]

  belongMat1 <- calcBelongMatrix(centers, data, 1.5)

  centers1 <- as.matrix(calcCentroids(data,belongMat1,1.5))
  centers2 <- as.matrix(calcSWFCCentroids(data, wdata, belongMat1, 1.5, 0))

  diff1 <- centers1 - centers2
  tot1 <- round(sum(diff1),8)
  expect_equal(tot1,0)
})
