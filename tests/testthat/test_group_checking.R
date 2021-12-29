context("group checking tests")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the euclidean distance function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The group matching function must returned what is expected",{

  mat1 <- rbind(
    c(0,1,0),
    c(0,1,0),
    c(1,0,0),
    c(1,0,0),
    c(0,0,1),
    c(0,0,1)
  )

  mat2 <- mat1[,c(3,1,2)]

  mat3 <- groups_matching(mat1, mat2)

  expect_equal(sum(abs(mat1-mat3)),0)

})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the clustering indices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Clustering indices should be able to distinguish good and bad clustering",{

  df1 <- data.frame(
    x = c(rnorm(50, mean = 0, sd = 0.2), rnorm(50, mean = 3, sd = 0.2)),
    y = c(rnorm(50, mean = 0, sd = 0.2), rnorm(50, mean = 3, sd = 0.2))
  )


  clus1 <- CMeans(df1, k = 2, m = 1.5)
  clus2 <- CMeans(df1, k = 5, m = 1.5)

  idxs_names <- c("FukuyamaSugeno.index","DaviesBoulin.index", "CalinskiHarabasz.index", "GD43.index", "GD53.index","Negentropy.index")
  upp_better <- c(
    FALSE, FALSE,TRUE,TRUE,TRUE,FALSE
  )

  idx1 <- calcqualityIndexes(data = df1,
                             belongmatrix = clus1$Belongings,
                             m = 1.5,
                             indices = idxs_names
                             )

  idx2 <- calcqualityIndexes(data = df1,
                             belongmatrix = clus2$Belongings,
                             m = 1.5,
                             indices = idxs_names
  )

  tests <- sapply(1:length(idx1), function(i){
    (idx1[[i]] > idx2[[i]]) == upp_better[[i]]
  })

  tot_test <- any(tests == FALSE)
  expect_false(tot_test)

})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the bootsraping function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The bootstraping function should be able to distinguish stable and unstable groups",{

  set.seed(123)
  df1 <- data.frame(
    x = c(rnorm(50, mean = 0, sd = 0.2), rnorm(50, mean = 3, sd = 0.2)),
    y = c(rnorm(50, mean = 0, sd = 0.2), rnorm(50, mean = 3, sd = 0.2))
  )


  clus1 <- CMeans(df1, k = 3, m = 1.5)
  # v2 should be the solid group
  bootValues <- boot_group_validation(clus1, nsim = 50)
  jacard_means <- colMeans(bootValues$group_consistency)

  tests <- c(jacard_means[[1]] < 0.6, jacard_means[[2]] > 0.7, jacard_means[[3]] < 0.6)
  expect_equal(sum(tests), 3)

})


test_that("The bootstraping function should be able to distinguish stable and unstable groups (multicore)",{

  set.seed(123)
  df1 <- data.frame(
    x = c(rnorm(50, mean = 0, sd = 0.2), rnorm(50, mean = 3, sd = 0.2)),
    y = c(rnorm(50, mean = 0, sd = 0.2), rnorm(50, mean = 3, sd = 0.2))
  )

  clus1 <- CMeans(df1, k = 3, m = 1.5)
  # v2 should be the solid group
  future::plan(future::multisession(workers=2))
  bootValues <- boot_group_validation.mc(clus1, nsim = 50, seed = 123)

  jacard_means <- colMeans(bootValues$group_consistency)

  tests <- c(jacard_means[[1]] < 0.6, jacard_means[[2]] > 0.7, jacard_means[[3]] < 0.6)
  expect_equal(sum(tests), 3)

})
