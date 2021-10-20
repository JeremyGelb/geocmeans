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
