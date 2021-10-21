context("c++ focal functions")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the focal euclidean function on a matrix
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the focal euclidean function on a matrix",{
  mat <- rbind(
    c(1,2,3),
    c(1,2,3),
    c(1,2,3)
  )

  W <- matrix(1,nrow = 3, ncol = 3)

  expected <- rbind(
    c(2,4,2),
    c(3,6,3),
    c(2,4,2)
  )

  obtained <-focal_euclidean_mat_window(mat, W)

  expect_equal(sum(abs(expected - obtained)),0)
})


test_that("Testing the focal euclidean function on a list of matrices",{
  mat <- rbind(
    c(1,2,3),
    c(1,2,3),
    c(1,2,3)
  )

  W <- matrix(1,nrow = 3, ncol = 3)

  list_mat <- list(mat,mat,mat)

  expected <- rbind(
    c(2,4,2),
    c(3,6,3),
    c(2,4,2)
  )

  expected <- expected*3

  obtained <-focal_euclidean_list(list_mat, W)

  expect_equal(sum(abs(expected - obtained)),0)

})
