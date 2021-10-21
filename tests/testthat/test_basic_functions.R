context("basic function tests")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the euclidean distance function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The euclidean distance function return what is expected",{
  mat <- matrix(4, ncol = 3, nrow = 3)
  vec <- c(2,6,3)
  dists <- calcEuclideanDistance(mat,vec)
  expected <- rep(sum((vec-4)**2), times = 3)
  expect_equal(dists, expected)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the cat_to_belongings
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("Testing the cat_to_belongings function",{
  values <- c(1,3,2,1)
  expected <- cbind(
    c(1,0,0,1),
    c(0,0,1,0),
    c(0,1,0,0)
  )
  obtained <- cat_to_belongings(values)
  expect_true(sum(obtained - expected) == 0)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the calcLaggedData
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_that("Testing the calcLaggedData function",{
  df1 <- data.frame(
    x = c(1,2,3,4)
  )

  nb_mat <- rbind(
    c(0,1,1,0),
    c(1,0,0,1),
    c(1,0,0,1),
    c(0,1,1,0)
  )

  nb <- spdep::mat2listw(nb_mat,style = "W")

  # for a simple mean
  obtained <- calcLaggedData(df1,nb)$x
  expected <- c(mean(c(2,3)), mean(c(1,4)), mean(c(1,4)), mean(c(2,3)))
  test1 <- (any(expected != obtained))==FALSE

  # for a median
  obtained <- calcLaggedData(df1,nb, method = "median")$x
  expected <- c(3,4,4,3)
  test2 <- (any(expected != obtained))==FALSE

  expect_true(test1 & test2)

})
