context("belonging matrices")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing that the obtained membership matrices have a rowSums = 1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Membership matrices must have a rowSums equal to 1",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")
  FCM <- CMeans(dataset, k = 5, m = 1.5, standardize = TRUE)
  GFCM <- GCMeans(dataset, k = 5, m = 1.5, beta = 0.5, standardize = TRUE)
  SFCM <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
  SGFCM <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = TRUE)
  x1 <- round(1-rowSums(FCM$Belongings),6)
  x2 <- round(1-rowSums(GFCM$Belongings),6)
  x3 <- round(1-rowSums(SFCM$Belongings),6)
  x4 <- round(1-rowSums(SGFCM$Belongings),6)
  expect_equal(sum(x1) + sum(x2) + sum(x3) + sum(x4),0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing that the algorithms have identical results if beta or alpha = 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("SFCM and FCM should yield identical results if alpha = 0 for SFCM",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")
  result1 <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 0, standardize = TRUE, seed = 123)
  result2 <- CMeans(dataset, k = 5, m = 1.5, standardize = TRUE, seed = 123)
  x <- round(result1$Belongings - result2$Belongings,8)
  expect_equal(sum(x),0)
})


test_that("GFCM and FCM should yield identical results if beta = 0 for GFCM",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]
  result1 <- GCMeans(dataset, k = 5, m = 1.5, beta = 0, standardize = TRUE, seed = 123)
  result2 <- CMeans(dataset, k = 5, m = 1.5, standardize = TRUE, seed = 123)
  x <- round(result1$Belongings - result2$Belongings,8)
  expect_equal(sum(x),0)
})


test_that("SGFCM and FCM should yield identical results if beta = 0 and alpha = 0 for SGFCM",{
  data(LyonIris)
  library(spdep)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]
  queen <- poly2nb(LyonIris,queen=TRUE)
  Wqueen <- nb2listw(queen,style="W")
  result1 <- SGFCMeans(dataset,Wqueen, k = 5, m = 1.5, alpha = 0, beta = 0, standardize = TRUE, seed = 123)
  result2 <- CMeans(dataset, k = 5, m = 1.5, standardize = TRUE, seed = 123)
  x <- round(result1$Belongings - result2$Belongings,8)
  expect_equal(sum(x),0)
})
