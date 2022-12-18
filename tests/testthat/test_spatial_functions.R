context("spatial function tests")
library(terra)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the global and local moran I
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The global moran I calculated is the same as the one in raster package",{
  set.seed(1254)
  mat <- matrix(sample(1:10, size = 100, replace = TRUE), nrow = 10, ncol = 10)
  w <- matrix(1, nrow = 3, ncol = 3)
  w[2,2] <- 0
  # raster::Moran(raster::raster(mat), w = w), expected value : 0.07608944
  expected <- 0.07608944
  val <- calc_moran_raster(terra::rast(mat), w = w)
  expect_equal(val, expected)
})


test_that("The local moran I calculated is the same as the one in raster package",{
  set.seed(1254)
  mat <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5, ncol = 5)
  #v1 <- raster::as.matrix(MoranLocal(raster::raster(mat)))
  expected <- c(-1.158590254, -0.008237818,  0.087379712, -0.005295740, -0.220802942)
  v2 <- terra::as.matrix(calc_local_moran_raster(terra::rast(mat), window = matrix(1, nrow = 3, ncol = 3)), wide = TRUE)
  expect_equal(expected, v2[1,])
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
  vals <- calcELSA(categories, nb, matdist = dists)
  testB <- round(vals[[5]],3) == 0.503

  ## test C
  categories <- rep(1, times = 9)
  categories[[3]] <- 2
  vals <- calcELSA(categories, nb, matdist = dists)
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
  vals <- calcFuzzyELSA(mat, nb, matdist = dists)
  testB <- round(vals[[5]],3) == 0.503

  ## test C
  mat <- matrix(0, nrow = 9, ncol =2)
  mat[,1] <- 1
  mat[3,] <- c(0,1)

  my_object <- FCMres(list(
    "Data" = data.frame(x = rep(0, times = 9),
                        y =rep(0, times = 9)),
    "Belongings" = mat,
    "Centers" = rbind(c(0,1),
                      c(0,0)),
    "m" = 1,
    "algo" = "cmeans"
  ))

  vals <- calcFuzzyELSA(my_object, nb, matdist = dists)
  testC <- round(vals[[5]],3) == 0.063

  ## test D
  dists <- rbind(c(0,1,1), c(1,0,1),c(1,1,0))
  mat <- cbind(
    c(1,0,1,0,0,0,1,0,1),
    c(0,1,0,1,0,1,0,1,0),
    c(0,0,0,0,1,0,0,0,0)
  )
  vals <- elsa_fuzzy_vector(mat, nb, matdist = dists)
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
  vals <- elsa_raster(mat, window, dists)
  testA <- round(vals[[5]],3) == 0

  ## test B
  mat <- matrix(0, ncol = 3, nrow = 3)
  mat[2,2] <- 1
  #vals <- raster::values(elsa_raster(raster::raster(mat), window, dists))
  vals <- terra::values(elsa_raster(terra::rast(mat),window, dists),mat = FALSE)
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
  vals <- terra::values(calcFuzzyELSA(object = list(
    terra::rast(mat1),
    terra::rast(mat2)),
                window = window, matdist = dists
  ))
  testB <- round(vals[[5]],3) == 0.503

  ## test C
  mat1 <- matrix(0, nrow = 3, ncol =3)
  mat2 <- matrix(1, nrow = 3, ncol =3)
  mat1[1,3] <- 1
  mat2[1,3] <- 0

  my_object <- FCMres(list(
    "Data" = list(terra::rast(mat1),terra::rast(mat1)),
    "rasters" = list(terra::rast(mat1),terra::rast(mat2)),
    "Centers" = rbind(c(0,1),
                      c(0,0)),
    "m" = 1,
    "algo" = "cmeans"
  ))

  vals <- terra::values(calcFuzzyELSA(object = my_object,
                        window = window, matdist = dists
  ), mat = FALSE)

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


  vals <- terra::values(calcFuzzyElsa_raster(list(
    terra::rast(mat1),
    terra::rast(mat2),
    terra::rast(mat3)),
    window = window, matdist = dists
  ), mat = FALSE)

  #arr <- array(c(mat1,mat2,mat3), c(3,3,3))
  #vals <- Elsa_fuzzy_matrix_window(arr, window, dists)
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing sp consistency index
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the spatial consistency index for a vector dataset",{

  # first, creating a membership matrix for 4 observations (squarred 2x2) and 2 groups
  belong_mat <- matrix(0,nrow = 4, ncol = 2)
  belong_mat[1,1] <- 1
  belong_mat[2,1] <- 1
  belong_mat[3,2] <- 1
  belong_mat[4,2] <- 1

  # second, creating the spatial neighbouring matrix
  neighmat <- matrix(0, ncol = 4, nrow = 4)
  neighmat[1,2] <- 1
  neighmat[2,1] <- 1
  neighmat[1,3] <- 1
  neighmat[3,1] <- 1

  neighmat[4,2] <- 1
  neighmat[2,4] <- 1
  neighmat[4,3] <- 1
  neighmat[3,4] <- 1
  nb <- spdep::mat2listw(neighmat)

  # each observation has a diff of 2 with a neighbour
  expected <- 4*2

  obtained <- spConsistency(belong_mat, nblistw = nb,nrep = 5)

  expect_equal(expected, obtained$sum_diff)

})


test_that("Testing the spatial consistency index for a raster dataset",{

  # first, creating a membership matrix for 4 observations (squarred 2x2) and 2 groups
  rast1 <- rbind(
    c(1,1),
    c(0,0)
  )
  rast1 <- terra::rast(rast1)

  rast2 <- rbind(
    c(0,0),
    c(1,1)
  )
  rast2 <- terra::rast(rast2)

  rasters <- list(rast1,rast2)
  W <- matrix(1, nrow = 3, ncol = 3)
  Data <- rasters
  centers <- data.frame(
    x1 = c(1,0),
    x2 = c(0,1)
  )

  myFCMres <- FCMres(list(
    "Data" = Data,
    "Centers" = centers,
    "rasters" = rasters,
    "m" = 1,
    "algo" = "kmeans"
  ))


  # each observation has a diff of 2 with two neighbours
  expected <- 2 * 2 * 4

  obtained <- spConsistency(myFCMres, window = matrix(1, nrow = 3, ncol = 3),nrep = 5)

  expect_equal(expected, obtained$sum_diff)

})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the adjusting of spatial weight by semantical distance
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Testing the function adjustSpatialWeights",{

  df1 <- data.frame(
    x = c(1,1,0,0),
    y = c(1,1,0,0)
  )

  neighmat <- rbind(
    c(0,1,1,0),
    c(1,0,0,1),
    c(1,0,0,1),
    c(0,1,1,0)
  )

  nb <- spdep::mat2listw(neighmat, style = "W")

  # expected weights
  distmat <- as.matrix(dist(df1)**2)
  distmat <- 1/distmat
  distmat[is.infinite(distmat)] <- 0
  distmat <- distmat * neighmat
  distmat <- distmat / rowSums(distmat)

  # obtained weight
  expect_warning({adj_nb <- adjustSpatialWeights(data = df1, listw = nb$neighbours, style = "W")})
  obtained <- round(spdep::listw2mat(adj_nb),5)
  diff <- sum(abs((distmat - obtained)))
  expect_equal(diff, 0)

})

