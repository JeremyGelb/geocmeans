context("raster functions")


test_that("Testing calcWdataRaster, it should return what is expected",{

  # preparing some test data
  w <- matrix(1, nrow = 3, ncol = 3)
  w[2,2] <- 0
  fun <- "mean"
  dataset <- matrix(rnorm(25,0,1), nrow = 5, ncol = 5)
  missing <- !is.na(c(dataset))

  obtained <- matrix(calcWdataRaster(w, list(terra::rast(dataset)), fun, missing), nrow = 5, ncol = 5)

  ob1 <- obtained[3,3]
  ex1 <- weighted.mean(dataset[2:4,2:4] , w)

  ob2 <- obtained[4,4]
  ex2 <-  weighted.mean(dataset[3:5,3:5] , w)

  ob3 <- obtained[1,1]
  ex3 <- weighted.mean(dataset[1:2,1:2] , rbind(c(0,1),c(1,1)))

  expect_equal(ob1, ex1)
  expect_equal(ob2, ex2)
  expect_equal(ob3, ex3)

})




test_that("Testing check_dim, it should raise an error when a window is not well specified",{

  w1 <- matrix(1, nrow = 3, ncol = 3)
  w2 <- matrix(1, nrow = 2, ncol = 2)
  w3 <- matrix(1, nrow = 3, ncol = 4)

  # this should not raise an error
  check_window(w1)

  # this should raise an error
  expect_error(check_window(w2))

  # this should also raise an error
  expect_error(check_window(w3))

})



test_that("Testing check_raters_dims, it should raise an error when a raster is not with the same dimensions",{

  w1 <- matrix(1, nrow = 3, ncol = 3)
  w2 <- matrix(1, nrow = 2, ncol = 2)
  w3 <- matrix(1, nrow = 3, ncol = 4)

  dataset <- list(terra::rast(w1), terra::rast(w2), terra::rast(w3))

  # this sould raise an error
  expect_error(check_raters_dims(dataset))

  # and this should not
  dataset <- list(terra::rast(w1), terra::rast(w1), terra::rast(w1))
  result <- check_raters_dims(dataset)

})
