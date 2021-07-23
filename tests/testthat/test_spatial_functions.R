context("spatial function tests")
library(raster)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the global and local moran I
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The global moran I calculated is the same as the one in raster package",{
  mat <- matrix(sample(1:10, size = 100, replace = TRUE), nrow = 10, ncol = 10)
  rast <- raster(mat)
  v1 <- Moran(rast)
  v2 <- moranI_matrix_window(mat, window = matrix(1, nrow = 3, ncol = 3))
  expect_equal(v1, v2)
})


test_that("The local moran I calculated is the same as the one in raster package",{
  mat <- matrix(sample(1:10, size = 100, replace = TRUE), nrow = 10, ncol = 10)
  rast <- raster(mat)
  v1 <- raster::values(MoranLocal(rast))
  v2 <- local_moranI_matrix_window(mat, window = matrix(1, nrow = 3, ncol = 3))
  expect_equal(v1, v2)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the ELSA function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### for vector

test_that("The ELSA index calculated on vector should give the same values as in the original paper",{

  categories <- rep(1, times = 9)
  categories[[5]] <- 2
  dists <- rbind(c(0,1), c(1,0))
  neighmat <- matrix(0, ncol = 9, nrow = 9)
  p1 <- c(1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,7,7,7,8,8,8,8,8,9,9,9)
  p2 <- c(2,4,5,1,4,5,6,3,2,5,6,1,2,5,7,8,1,2,3,4,6,7,8,9,2,3,5,8,9,4,5,8,4,5,6,7,9,5,6,8)

  for(i in 1:length(p1)){
    neighmat[p1[[i]], p2[[i]]] <- 1
  }
  nb <- spdep::mat2listw(neighmat)

  ## testA
  vals <- elsa_vector(rep(1,9), nb, dists)
  testA <- round(vals[[5]],3) == 0

  ## test B
  vals <- elsa_vector(categories, nb, dists)
  testB <- round(vals[[5]],3) == 0.503

  ## test C
  categories <- rep(1, times = 9)
  categories[[3]] <- 2
  vals <- elsa_vector(categories, nb, dists)
  testC <- round(vals[[5]],3) == 0.063

  ## test D
  dists <- rbind(c(0,1,1), c(1,0,1),c(1,1,0))
  categories <- c(1,2,1,2,3,2,1,2,1)
  vals <- elsa_vector(categories, nb, dists)
  testD <- round(vals[[5]],3) == 0.878

  ## test E
  dists <- rbind(c(0,1,1), c(1,0,1),c(1,1,0))
  categories <- c(2,2,2,1,3,1,2,2,2)
  vals <- elsa_vector(categories, nb, dists)
  testE <- round(vals[[5]],3) == 0.773

  expect_true(testA & testB & testC & testD & testE)
})


test_that("The fuzzy ELSA index calculated on vector should give the same values as in the original paper",{

  dists <- rbind(c(0,1), c(1,0))
  neighmat <- matrix(0, ncol = 9, nrow = 9)
  p1 <- c(1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,7,7,7,8,8,8,8,8,9,9,9)
  p2 <- c(2,4,5,1,4,5,6,3,2,5,6,1,2,5,7,8,1,2,3,4,6,7,8,9,2,3,5,8,9,4,5,8,4,5,6,7,9,5,6,8)

  for(i in 1:length(p1)){
    neighmat[p1[[i]], p2[[i]]] <- 1
  }
  nb <- spdep::mat2listw(neighmat)

  ## testA
  mat <- matrix(0, nrow = 9, ncol =2)
  mat[,1] <- 1
  vals <- elsa_fuzzy_vector(mat, nb, dists)
  testA <- round(vals[[5]],3) == 0

  ## test B
  mat <- matrix(0, nrow = 9, ncol =2)
  mat[,1] <- 1
  mat[5,] <- c(0,1)
  vals <- elsa_fuzzy_vector(mat, nb, dists)
  testB <- round(vals[[5]],3) == 0.503

  ## test C
  mat <- matrix(0, nrow = 9, ncol =2)
  mat[,1] <- 1
  mat[3,] <- c(0,1)
  vals <- elsa_fuzzy_vector(mat, nb, dists)
  testC <- round(vals[[5]],3) == 0.063

  ## test D
  dists <- rbind(c(0,1,1), c(1,0,1),c(1,1,0))
  mat <- cbind(
    c(1,0,1,0,0,0,1,0,1),
    c(0,1,0,1,0,1,0,1,0),
    c(0,0,0,0,1,0,0,0,0)
  )
  vals <- elsa_fuzzy_vector(mat, nb, dists)
  testD <- round(vals[[5]],3) == 0.878

  ## test E
  dists <- rbind(c(0,1,1), c(1,0,1),c(1,1,0))
  mat <- cbind(
    c(1,1,1,0,0,0,1,1,1),
    c(0,0,0,1,0,1,0,0,0),
    c(0,0,0,0,1,0,0,0,0)
  )
  vals <- elsa_fuzzy_vector(mat, nb, dists)
  testE <- round(vals[[5]],3) == 0.773

  expect_true(testA & testB & testC & testD & testE)
})


### for raster

test_that("The ELSA index calculated on raster should give the same values as in the original paper",{

  dists <- rbind(c(0,1), c(1,0))
  window <- matrix(1, ncol = 3, nrow = 3)

  ## testA
  mat <- matrix(0, ncol = 3, nrow = 3)
  vals <- Elsa_categorical_matrix_window(mat, window, dists)
  testA <- round(vals[[5]],3) == 0

  ## test B
  mat <- matrix(0, ncol = 3, nrow = 3)
  mat[2,2] <- 1
  vals <- Elsa_categorical_matrix_window(mat, window, dists)
  testB <- round(vals[[5]],3) == 0.503

  ## test C
  mat <- matrix(0, ncol = 3, nrow = 3)
  mat[1,3] <- 1
  vals <- Elsa_categorical_matrix_window(mat, window, dists)
  testC <- round(vals[[5]],3) == 0.063

  ## test D
  dists <- rbind(c(0,1,1), c(1,0,1),c(1,1,0))
  mat <- matrix(c(0,1,0,1,2,1,0,1,0), ncol = 3, nrow = 3)
  vals <- Elsa_categorical_matrix_window(mat, window, dists)
  testD <- round(vals[[5]],3) == 0.878

  ## test E
  dists <- rbind(c(0,1,1), c(1,0,1),c(1,1,0))
  mat <- matrix(c(0,0,0,2,1,2,0,0,0), ncol = 3, nrow = 3)
  vals <- Elsa_categorical_matrix_window(mat, window, dists)
  testE <- round(vals[[5]],3) == 0.773

  expect_true(testA & testB & testC & testD & testE)
})


test_that("The fuzzy ELSA index calculated on raster should give the same values as in the original paper",{

  dists <- rbind(c(0,1), c(1,0))
  mat <- matrix(0, ncol = 3, nrow = 3)
  window <- matrix(1, ncol = 3, nrow = 3)

  ## testA
  mat1 <- matrix(0, nrow = 3, ncol =3)
  mat2 <- matrix(1, nrow = 3, ncol =3)
  arr <- array(c(mat1,mat2), c(3,3,2))
  vals <- Elsa_fuzzy_matrix_window(arr, window, dists)
  testA <- round(vals[[5]],3) == 0

  ## test B
  mat1 <- matrix(0, nrow = 3, ncol =3)
  mat2 <- matrix(1, nrow = 3, ncol =3)
  mat1[2,2] <- 1
  mat2[2,2] <- 0
  arr <- array(c(mat1,mat2), c(3,3,2))
  vals <- Elsa_fuzzy_matrix_window(arr, window, dists)
  testB <- round(vals[[5]],3) == 0.503

  ## test C
  mat1 <- matrix(0, nrow = 3, ncol =3)
  mat2 <- matrix(1, nrow = 3, ncol =3)
  mat1[1,3] <- 1
  mat2[1,3] <- 0
  arr <- array(c(mat1,mat2), c(3,3,2))
  vals <- Elsa_fuzzy_matrix_window(arr, window, dists)
  testC <- round(vals[[5]],3) == 0.063

  ## test D
  dists <- matrix(1, nrow = 3, ncol = 3)
  diag(dists) <- 0
  mat1 <- cbind(c(0,0,0),
                c(0,1,0),
                c(0,0,0)
                )
  mat2 <- cbind(c(0,1,0),
                c(1,0,1),
                c(0,1,0)
  )
  mat3 <- cbind(c(1,0,1),
                c(0,0,0),
                c(1,0,1)
  )
  arr <- array(c(mat1,mat2,mat3), c(3,3,3))
  vals <- Elsa_fuzzy_matrix_window(arr, window, dists)
  testD <- round(vals[[5]],3) == 0.878

  ## test E
  mat1 <- cbind(c(0,0,0),
                c(0,1,0),
                c(0,0,0)
  )
  mat2 <- cbind(c(1,1,1),
                c(0,0,0),
                c(1,1,1)
  )
  mat3 <- cbind(c(0,0,0),
                c(1,0,1),
                c(0,0,0)
  )
  arr <- array(c(mat1,mat2, mat3), c(3,3,3))
  vals <- Elsa_fuzzy_matrix_window(arr, window, dists)
  testE <- round(vals[[5]],3) == 0.773

  expect_true(testA & testB & testC & testD & testE)
})

