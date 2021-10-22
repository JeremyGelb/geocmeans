#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the function to evalaute FCM parameters
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The evaluate parameters function must return the same values as calculated by the simple functions",{

  ## situation
  data(LyonIris)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
  "TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]

  ks <- c(2,3)
  ms <- c(1.5,2)
  myseed <- 123

  idxs <- c("Silhouette.index", "Partition.entropy", "Partition.coeff", "XieBeni.index", "Explained.inertia")
  list_res <- list()
  for (k in ks){
    for (m in ms){
      res <- CMeans(dataset, k, m, seed = myseed, standardize = FALSE, verbose =F, tol = 0.001)
      vals <- calcqualityIndexes(dataset,res$Belongings,m = m,indices = idxs)
      list_res[[length(list_res)+1]] <- unlist(vals)
    }
  }
  df1 <- data.frame(do.call(rbind, list_res))
  df1$k <- c(2,2,3,3)
  df1$m <- c(1.5,2,1.5,2)
  df1$oid <- paste(df1$k,df1$m, sep = "_")
  df1 <- df1[order(df1$oid),]

  values <- select_parameters("FCM", dataset, k = c(2,3), m = c(1.5,2),
      spconsist=FALSE, seed = 123, standardize = F, tol = 0.001, indices = idxs)
  values$oid <- paste(values$k, values$m, sep = "_")
  values <- values[order(values$oid),]

  obtained <- sum(round(values[,1:5],3) - round(df1[,1:5],3))
  expect_equal(obtained, 0)
})


test_that("The evaluate parameters function (multicore) must return the same values as calculated by the simple functions",{

  ## situation
  data(LyonIris)
  AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
                     "TxChom1564","Pct_brevet","NivVieMed")
  dataset <- LyonIris@data[AnalysisFields]

  ks <- c(2,3)
  ms <- c(1.5,2)
  myseed <- 123

  idxs <- c("Silhouette.index", "Partition.entropy", "Partition.coeff", "XieBeni.index", "Explained.inertia")
  list_res <- list()
  for (k in ks){
    for (m in ms){
      res <- CMeans(dataset, k, m, seed = myseed, standardize = FALSE, verbose =F, tol = 0.000001)
      vals <- calcqualityIndexes(dataset,res$Belongings,m = m,indices = idxs)
      list_res[[length(list_res)+1]] <- unlist(vals)
    }
  }
  df1 <- data.frame(do.call(rbind, list_res))
  df1$k <- c(2,2,3,3)
  df1$m <- c(1.5,2,1.5,2)
  df1$oid <- paste(df1$k,df1$m, sep = "_")
  df1 <- df1[order(df1$oid),]

  future::plan(future::multiprocess(workers=2))
  values <- select_parameters.mc("FCM", dataset, k = c(2,3), m = c(1.5,2),
                              spconsist=FALSE, seed = 123, standardize = F, tol = 0.000001,
                              indices = idxs)
  values$oid <- paste(values$k, values$m, sep = "_")
  values <- values[order(values$oid),]

  obtained <- sum(round(values[,1:5],3) - round(df1[,1:5],3))
  expect_equal(obtained, 0)
})
