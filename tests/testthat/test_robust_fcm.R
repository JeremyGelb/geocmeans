context("robust and noisy FCM tests")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the function producing the sigmas
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("The calcRobustSigmas function must returned what is expected",{

 data <- rbind(
   c(1,2,3),
   c(1,1,1),
   c(1,2,2),
   c(1,1,1),
   c(1,2,3)
 )
 belong_matrix <- rbind(
   c(0.2,0.8),
   c(0.4,0.6),
   c(0.5,0.5),
   c(0.5,0.5),
   c(0.7,0.3)
 )
 m <- 1.5

 centers <- calcCentroids(data, belong_matrix, m)
 powered <- belong_matrix**m

 expected <- sqrt(sapply(1:2, function(i){
   ui <- powered[,i]
   denom <- sum(ui)
   num <- sum(calcEuclideanDistance2(data, centers[i,]) * ui)
   return(num/denom)
 }))

 obtained <- calcRobustSigmas(data, belong_matrix, centers, m)

 expect_equal(expected, obtained)

})
