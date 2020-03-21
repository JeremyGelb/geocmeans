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
