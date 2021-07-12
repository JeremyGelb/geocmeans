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
